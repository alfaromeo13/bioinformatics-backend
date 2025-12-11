#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re, io, sys, logging, logging.handlers, argparse
from pathlib import Path
from typing import List, Optional, Dict
import multiprocessing as mp

# === your modules ===
from scripts.residue_detection import ResidueDetector
from scripts.ligand_prep import LigandPreparation
from scripts.prepare_pdb import PreparePdb
from scripts.step1 import Step1
from scripts.step2 import Step2
from scripts.step3_inte import Step3INTE
from scripts.step3_gbmv import Step3GBMV
from scripts.heatmap import HeatmapGenerator
from scripts.cleanup import Cleanup
from scripts.error_handle import Errors
from scripts.greeting import Greeting

# ---------------- logging (file + console, capture prints) ----------------
def init_logging(logfile: Optional[str], verbosity: int = 1):
    level = logging.WARNING if verbosity == 0 else logging.INFO if verbosity == 1 else logging.DEBUG
    handlers = [logging.StreamHandler(sys.stdout)]
    if logfile:
        Path(logfile).parent.mkdir(parents=True, exist_ok=True)
        handlers.append(logging.handlers.RotatingFileHandler(
            logfile, maxBytes=5_000_000, backupCount=3, encoding="utf-8"
        ))
    logging.basicConfig(
        level=level,
        format="%(asctime)s | %(levelname)-8s | %(name)s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        handlers=handlers,
        force=True,
    )
    logging.captureWarnings(True)
    def _exhook(exc_type, exc, tb):
        logging.getLogger("FATAL").exception("Uncaught exception", exc_info=(exc_type, exc, tb))
        sys.__excepthook__(exc_type, exc, tb)
    sys.excepthook = _exhook

    class _Tee(io.TextIOBase):
        def __init__(self, stream, logger, level):
            self._stream = stream; self._logger = logger; self._level = level
        def write(self, buf):
            self._stream.write(buf)
            for line in buf.rstrip().splitlines():
                if line.strip():
                    self._logger.log(self._level, line)
            return len(buf)
        def flush(self): self._stream.flush()
    sys.stdout = _Tee(sys.stdout, logging.getLogger("STDOUT"), logging.INFO)
    sys.stderr = _Tee(sys.stderr, logging.getLogger("STDERR"), logging.ERROR)


def log_run_header(args):
    """Log basic system info and all CLI input parameters."""
    import datetime, platform
    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    logging.info("=" * 70)
    logging.info(f"Run started: {now}")
    logging.info(f"Working directory: {Path.cwd()}")
    logging.info(f"Python: {platform.python_version()} ({sys.executable})")
    logging.info(f"System: {platform.system()} {platform.release()} ({platform.machine()})")
    logging.info(f"CPU count: {mp.cpu_count()}")
    logging.info("Input parameters:")
    for key, value in vars(args).items():
        logging.info(f"  {key:20s}: {value}")
    logging.info("=" * 70)

# ---------------- helpers ----------------
_MUT_RE = re.compile(r"^(?P<wtres>[A-Z]{3})(?P<chain>[A-Za-z])(?P<resn>\d+)$")

def ensure_exists(path: Path, what: str):
    if not path.exists():
        raise FileNotFoundError(f"{what} not found: {path}")

def _normalize_mut(tok: str) -> Optional[str]:
    t = tok.replace(",", "").strip().upper()
    m = _MUT_RE.match(t)
    if not m:
        return None
    wt, ch, rn = m.group("wtres"), m.group("chain"), m.group("resn")
    if wt == "HIS":  # CHARMM default tautomer
        wt = "HSD"
    return f"{wt}{ch}{rn}"

def parse_mutations_arg(arg: Optional[str]) -> Optional[List[str]]:
    if not arg:
        return None
    toks = re.split(r"[,\s]+", arg.strip())
    out, seen = [], set()
    for tok in toks:
        if not tok:
            continue
        nm = _normalize_mut(tok)
        if not nm:
            logging.warning(f"Skipping malformed mutation token: {tok!r} (expected e.g. ALAC58)")
            continue
        if nm in seen:
            continue
        seen.add(nm); out.append(nm)
    return out or None

def filter_terminals_in_memory(pdb_path: str, mutations: List[str]) -> List[str]:
    """
    Remove mutations that target first/last residue of their chain (no setup file used).
    Lightweight PDB ATOM/HETATM parser; tolerates insertion codes by using numeric part.
    """
    per_chain = {}  # chain -> set of residue numbers
    with open(pdb_path, "r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            if len(line) < 26:
                continue
            chain = line[21].strip() or " "
            resseq = line[22:26].strip()
            m = re.match(r"(\d+)", resseq)
            if not m:
                continue
            resn = int(m.group(1))
            per_chain.setdefault(chain, set()).add(resn)
    chain_minmax = {ch: (min(v), max(v)) for ch, v in per_chain.items() if v}

    kept, dropped = [], []
    for mut in mutations:
        ch, rn = mut[3], int(mut[4:])
        mnmx = chain_minmax.get(ch)
        if mnmx and (rn == mnmx[0] or rn == mnmx[1]):
            dropped.append(mut)
        else:
            kept.append(mut)
    if dropped:
        logging.info(f"Removed terminal residues (in-memory): {', '.join(dropped)}")
    logging.info(f"{len(kept)} mutations remain after terminal filtering (in-memory).")
    return kept

def residue_detection_cli_only(pdb_name: str, cutoff: float,
                               protein_chains: List[str], partner_chains: List[str]) -> List[str]:
    """
    Runs ResidueDetector WITHOUT a setup file.
    Requires ResidueDetector to expose 'detected_mutations' (list of tokens like ALAC58).
    """
    logging.info("Auto-detecting contact residues (CLI-only)…")
    rd = ResidueDetector(pdb_name, cutoff, protein_chains, partner_chains, None)
    rd.detect_distances()
    rd.check_mutations()
    muts = getattr(rd, "detected_mutations", None)
    if not muts:
        logging.error("ResidueDetector must expose 'detected_mutations' when used without a setup file.")
        raise RuntimeError("Automatic detection failed (no 'detected_mutations'). Provide --mutations explicitly.")
    norm = [_normalize_mut(m) for m in muts]
    return [m for m in norm if m]

# ---------------- pipeline wrappers ----------------
def do_ligprep(cgenff: Optional[str]):
    if cgenff:
        lp = LigandPreparation(cgenff)
        lp.split_cgen_out(); lp.hydrogens_under_carbons()
        logging.info(f"CGENFF custom parameters adjusted: {cgenff}")
    else:
        logging.info("No CGENFF file (skipping ligand preparation).")

def do_prepare_protein(pdb_name: str):
    pp = PreparePdb(pdb_name)
    pp.hetatm2atm(); pp.his2hsd(); pp.replace_ion(); pp.replace_dna_rna()
    pp.chain_spliter(); pp.small_caps()

def do_step1(pdb_name: str, charmm_dir: str):
    s1 = Step1(pdb_name, charmm_dir)
    s1.step1_maker(); s1.run_step1()

def do_step2(mutation_list: List[str], charmm_dir: str, thread_num: int):
    for mut in mutation_list:
        wtres, ch, resn = mut[:3], mut[3], int(mut[4:])
        s2 = Step2(ch, 'PRO'+ch, resn, wtres, charmm_dir, thread_num)
        s2.generate_charmm_input(); s2.run_step2()

def do_step3_inte(mutation_list: List[str], params: Dict):
    for mut in mutation_list:
        wtres, ch, resn = mut[:3], mut[3], int(mut[4:])
        s3 = Step3INTE(params["cgenff"], params["abnr"], 'PRO'+ch, resn, wtres,
                       params["protein_chains"], params["partner_chains"],
                       params["charmm_dir"], params["pdb"], params["threads"])
        s3.generate_all_chains(); s3.generate_non_mut_chid()
        s3.create_step3_inte_input(); s3.run_step3()

def do_step3_gbmv(mutation_list: List[str], params: Dict):
    for mut in mutation_list:
        wtres, ch, resn = mut[:3], mut[3], int(mut[4:])
        s3 = Step3GBMV(params["cgenff"], params["abnr"], 'PRO'+ch, resn, wtres,
                       params["protein_chains"], params["partner_chains"],
                       params["charmm_dir"], params["pdb"], params["threads"])
        s3.generate_all_chains(); s3.generate_non_mut_chid()
        s3.create_step3_gbmv_comp(); s3.run_step3_gbmv_comp()
    for mut in mutation_list:
        wtres, ch, resn = mut[:3], mut[3], int(mut[4:])
        s3 = Step3GBMV(params["cgenff"], params["abnr"], 'PRO'+ch, resn, wtres,
                       params["protein_chains"], params["partner_chains"],
                       params["charmm_dir"], params["pdb"], params["threads"])
        s3.calculate_gbmv()

def do_heatmap(mutation_list: List[str], pdb: str, method_label: str):
    aa = ['A','R','N','D','C','E','Q','G','H','I','L','K','M','F','P','S','T','W','Y','V']
    matrix_norm, matrix_abs = [], []
    for mut in mutation_list:
        ch, resn, wtres = mut[3], int(mut[4:]), mut[:3]
        hg = HeatmapGenerator(ch, resn, wtres, method_label)
        hg.inter_ener_data(); hg.get_ref_ener_lines()
        diff, diff_norm = hg.create_ener_diff_dict()
        matrix_abs.append([diff.get(a,0) for a in aa])
        matrix_norm.append([diff_norm.get(a,0) for a in aa])
    HeatmapGenerator(ch, resn, wtres, method_label).make_heatmap(matrix_norm, matrix_abs, mutation_list, aa, pdb, method_label)

def do_prechecks(protein_chains: List[str], partner_chains: List[str]):
    e = Errors(protein_chains, partner_chains)
    e.generate_all_chains()
    e.check_inte_crd()

def do_cleanup(pdb: str):
    c = Cleanup(pdb)
    c.zip_files(); c.delete_files()

# ---------------- CLI ----------------
def build_cli() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="In silico mutagenesis (CLI-only)")
    # required experiment inputs
    p.add_argument("--charmm-dir", required=True, type=str, help="Path to CHARMM executable")
    p.add_argument("--pdb", required=True, help="PDB file (e.g., 1brs.pdb)")
    p.add_argument("--protein-chains", required=True, help="Chains to mutate, e.g., 'C' or 'A,B'")
    p.add_argument("--partner-chains", required=True, help="Interacting partner chains, e.g., 'F' or 'J,L'")
    # mutagenesis controls
    p.add_argument("--mutations", help="Space/comma-separated list like 'ARGC59,GLUC60'; omit to auto-detect by cutoff")
    p.add_argument("--cutoff", type=float, default=5.0, help="Contact cutoff (Å) for auto-detection (default: 5.0)")
    p.add_argument("--abnr", type=int, default=50, help="ABNR minimization steps (default: 50)")
    p.add_argument("--cgenff", help="CGENFF .str/.rtf/.prm bundle for small molecules (optional)")
    p.add_argument("--threads", type=int, default=mp.cpu_count(), help="Thread count for Step2/3 (default: all CPUs)")
    # run control
    p.add_argument("--method", choices=["inte", "gbmv", "both"], default="both", help="Energy method(s) to run")
    p.add_argument("--skip-ligprep", action="store_true", help="Skip ligand preparation even if --cgenff is provided")
    p.add_argument("--skip-heatmap", action="store_true", help="Skip heatmap generation")
    p.add_argument("--skip-cleanup", action="store_true", help="Skip zip/delete outputs at the end")
    p.add_argument("--start-at", choices=["prep", "step1", "step2", "step3", "heatmap"], default="prep", help="Resume from stage")
    p.add_argument("--stop-after", choices=["prep", "step1", "step2", "step3", "heatmap", "done"], default="done", help="Stop after stage")
    # logging
    p.add_argument("-l", "--log", dest="logfile", help="Path to logfile (rotates at 5MB, keeps 3 backups)")
    p.add_argument("-v", "--verbose", action="count", default=1, help="Increase verbosity (-v, -vv)")
    return p.parse_args()

# ---------------- main ----------------
def main():
    args = build_cli()

    Greeting().start()

    init_logging(args.logfile, args.verbose)
    log_run_header(args)

    pdb_path = Path(args.pdb); ensure_exists(pdb_path, "PDB file")
    charmm_dir = args.charmm_dir  # CLI-only
    protein_chains = [c.strip() for c in re.split(r"[,\s]+", args.protein_chains) if c.strip()]
    partner_chains = [c.strip() for c in re.split(r"[,\s]+", args.partner_chains) if c.strip()]
    if not protein_chains or not partner_chains:
        raise ValueError("Both --protein-chains and --partner-chains must be non-empty.")



    # ===== PREP =====
    if args.start_at == "prep":
        # a) mutations: CLI list or b) auto-detect
        muts = parse_mutations_arg(args.mutations)
        if not muts:
            logging.info("No --mutations provided → automatic residue detection will be used.")
            muts = residue_detection_cli_only(str(pdb_path), args.cutoff, protein_chains, partner_chains)

        # terminal filtering in-memory (no setup reads/writes)
        muts = filter_terminals_in_memory(str(pdb_path), muts)
        if not muts:
            logging.error("No mutations remain after terminal filtering. Exiting.")
            sys.exit(2)

        if not args.skip_ligprep:
            do_ligprep(args.cgenff)
        else:
            logging.info("LigandPreparation skipped by flag.")

        do_prepare_protein(str(pdb_path))
        if args.stop_after == "prep":
            return

    # carry forward
    final_muts = muts
    params = dict(
        cgenff=args.cgenff,
        abnr=str(args.abnr),
        protein_chains=protein_chains,
        partner_chains=partner_chains,
        charmm_dir=charmm_dir,
        pdb=str(pdb_path),
        threads=int(args.threads),
    )

    # ===== STEP1 =====
    if args.start_at in ("prep", "step1"):
        do_step1(str(pdb_path), charmm_dir)
        if args.stop_after == "step1":
            return

    # ===== STEP2 =====
    if args.start_at in ("prep", "step1", "step2"):
        do_step2(final_muts, charmm_dir, int(args.threads))
        if args.stop_after == "step2":
            return

    # ===== PRECHECKS =====
    do_prechecks(protein_chains, partner_chains)

    # ===== STEP3 =====
    if args.start_at in ("prep", "step1", "step2", "step3"):
        if args.method in ("inte", "both"):
            do_step3_inte(final_muts, params)
        if args.method in ("gbmv", "both"):
            do_step3_gbmv(final_muts, params)
        if args.stop_after == "step3":
            return

    # ===== HEATMAP =====
    if not args.skip_heatmap and args.start_at in ("prep", "step1", "step2", "step3", "heatmap"):
        label = "gbmv" if args.method in ("gbmv", "both") else "inter"
        do_heatmap(final_muts, str(pdb_path), label)

    # ===== CLEANUP =====
    if not args.skip_cleanup and args.stop_after == "done":
        do_cleanup(str(pdb_path))

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        logging.exception(f"Fatal error: {e}")
        sys.exit(1)
