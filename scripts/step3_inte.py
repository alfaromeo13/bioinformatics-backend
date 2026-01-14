import os
import subprocess
import sys
import time
import threading
from collections import defaultdict


class Step3INTE:
    def __init__(self, cgenff_name, complex_abnr, pro_mut_chid, resn, wtres,protein_chains_list, ligand_chains_list,charmm_dir, pdb_name, thread_num):                 
        self.cgenff_name = cgenff_name
        self.pdb_name = pdb_name
        self.complex_abnr = complex_abnr
        self.pro_mut_chid = pro_mut_chid
        self.resn = resn
        self.wtres = wtres
        self.protein_chains_list = protein_chains_list
        self.ligand_chains_list = ligand_chains_list
        self.thread_num = thread_num
        self.A1 = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
                   'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HSD': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M',
                   'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
                   'TYR': 'Y', 'VAL': 'V'}
        self.all_chains = []
        self.non_mut_chid = []
        self.charmm_dir = charmm_dir
        self.subprocess_failed = False

    def disu_patch(self, wt, value):
        print("!Patch disulfide bonds based on SSBOND info in PDB")
        checker = self.pro_mut_chid + " " + str(self.resn)
        with open(self.pdb_name, 'r') as file:
            lines = file.readlines()

        patch_disu_lines = [] 

        for line in lines:
            if line.startswith('SSBOND'):
                fields = line.split()
                res1_chain = fields[3]
                res1_num = fields[4]
                res2_chain = fields[6]
                res2_num = fields[7]

                patch_disu_line = f"patch disu PRO{res1_chain} {res1_num} PRO{res2_chain} {res2_num} setup warn"
                patch_disu_lines.append(patch_disu_line)

        if any(checker in line for line in patch_disu_lines):
            if wt.lower() != 'c' or value.lower() != 'c':  # Convert wt and value to lowercase
                patch_disu_lines = [line for line in patch_disu_lines if checker not in line]

        for line in patch_disu_lines:
            print(line)
        #tega le če so SS bond, zaenkrat  ne vpliva        
        #print("autogenerate angl dihe\n")

    def generate_all_chains(self):
        for chain in self.protein_chains_list:
            self.all_chains.append('PRO' + chain)
        for chain in self.ligand_chains_list:
            self.all_chains.append('PRO' + chain)

    def generate_non_mut_chid(self):
        mut_chid = self.pro_mut_chid
        self.non_mut_chid = [chid for chid in self.all_chains if chid != mut_chid]
        for i in range(len(self.non_mut_chid)):
            self.non_mut_chid[i] = 'pro' + self.non_mut_chid[i]

    def process_non_mut(self):
        for chain_name in self.non_mut_chid:
            non_mut_chain_name = chain_name[3:]

            print("read psf card append name {}.psf".format(non_mut_chain_name.lower()))            
            print("read coor card append name {}.crd".format(non_mut_chain_name.lower()))              
    
    def charmm_inte(self):
        prot_chains = self.protein_chains_list
        lig_chains = self.ligand_chains_list

        output = 'inte sele segid pro{} '.format(prot_chains[0].lower())
        for chain in prot_chains[1:]:
            output += '.or. segid pro{} '.format(chain.lower())
        output += 'show end sele segid pro{} '.format(lig_chains[0].lower())
        for chain in lig_chains[1:]:
            output += '.or. segid pro{} '.format(chain.lower())
        output += 'show end'
        print(output)

    #ta funkcija za inte izbere mutacijo ter druge chaine ki so navedeni kot mutprot
    def charmm_iso_inte_mutprot(self):
        prot_chains = self.protein_chains_list
        lig_chains = self.ligand_chains_list
        pro_mut_chid = self.pro_mut_chid[3]
        mut_residue = self.resn

        output = 'inte sele '
        prot_processed = False
        lig_processed = False

        for i, chain in enumerate(prot_chains):
            if i > 0:
                output += '.or. '
            
            if chain == pro_mut_chid:
                output += 'segid pro{} .and. resi {} '.format(chain.lower(), mut_residue)
            else:
                output += 'segid pro{} '.format(chain.lower())
            
            if chain == prot_chains[-1] and not lig_processed:
                prot_processed = True
                output += 'end '
                output += 'sele '

        for i, chain in enumerate(lig_chains):
            if i > 0 and (prot_processed or i > 1):
                output += '.or. '
            
            if chain == pro_mut_chid:
                output += 'segid pro{} .and. resi {} '.format(chain.lower(), mut_residue)
            else:
                output += 'segid pro{} '.format(chain.lower())
            
            if i == len(lig_chains) - 1:
                lig_processed = True
                if prot_processed:
                    output += 'end '
        
        print(output)

    def ligand_rtf_prm(self):
        if self.cgenff_name is not None:
            print("open read card unit 10 name ligand.rtf")
            print("read rtf card unit 10 append")
            print("open read card unit 20 name ligand.prm")
            print("read para flex card unit 20 append")

    def charmm_join_inter_ener(self, wt):
        for value in self.A1.values():
            self.process_non_mut()
            print("read psf card append name {}_{}_{}2{}.psf".format(self.pro_mut_chid, self.resn, wt, value))
            print("read coor card append name {}_{}_{}2{}.crd".format(self.pro_mut_chid, self.resn, wt, value))
            #patchanje disulfidnih vezi
            self.disu_patch(wt, value)
            # MINIMIZATION
            #nonbondi, ki so isti kot pri uporabi GBMV!
            print("ENERGY atom CDIE eps 1 cutnb 21 ctofnb 18 ctonnb 16 switch vswitch\n")
            print("mini SD nstep 50 nprint 50")
            print("mini ABNR nstep {} nprint 50".format(self.complex_abnr))
            print("define msit  sele segid {} .and. resid {} end".format(self.pro_mut_chid,self.resn))
            print("cons fix sele .not. ( msit  .around. 10.0 ) end")
            print("mini ABNR nstep 10000 tolgrd 0.05 nprint 50\n")
            print("\ncons fix sele none end")
            print("ENERGY atom CDIE eps 1 cutnb 21 ctofnb 18 ctonnb 16 switch vswitch\n")
            # writing joined crd psf and pdb
            print("open write unit 10 card name joined_{}_{}_{}2{}.psf".format(self.pro_mut_chid, self.resn, wt, value))
            print("write psf unit 10 card")
            print("open write card unit 10 card name joined_{}_{}_{}2{}.crd".format(self.pro_mut_chid, self.resn, wt, value))
            print("write coor unit 10 card")
            print("open write card unit 10 card name joined_{}_{}_{}2{}.pdb".format(self.pro_mut_chid, self.resn, wt, value))
            print("write coor pdb unit 10 card")
            #Tole ne rabimo če gbmv dela
            print("set t {}_{}_{}2{}\n".format(self.pro_mut_chid, self.resn, wt, value))
            #Nonbond nerabimo, ker default dobimo iz ABNR, tam je tudi update komanda ki je inte nima avtomatsko
            self.charmm_inte()
            print("\nset e1 ?ener")
            mut_chid = self.pro_mut_chid[3:]
            print("open write unit 21 form name inter_ener_{}{}{}.dat append".format(self.wtres, self.resn, mut_chid).lower())
            print("write title unit 21")
            print("*@t          @e1")
            print("*\n")
            print("close unit 21")
            print("delete atom select all end")


    def create_step3_inte_input(self):
        mut_chid = self.pro_mut_chid[3:]
        with open('step3_{}{}{}_inte.inp'.format(self.wtres, self.resn, mut_chid).lower(), "w") as f:
            sys.stdout = f

            print("*GENERATED BY SEBALAB")
            print("*SMEHeC: Swift Mutational Energy Heatmap Calculator")
            print("*Joining and minimization of chains with mutations")
            print("*")
            print("stream charmm_dat/toppar.str")
            self.ligand_rtf_prm()
            mut_chid = self.pro_mut_chid[3:]
            print("open write unit 21 form name inter_ener_{}{}{}.dat".format(self.wtres, self.resn, mut_chid).lower())
            print("write title unit 21\n")
            
            wt = self.A1.get(self.wtres)
            self.charmm_join_inter_ener(wt)

            print("stop")
            sys.stdout = sys.__stdout__

    def run_step3(self):
        mut_chid = self.pro_mut_chid[3:].lower()
        mut_chid_caps = self.pro_mut_chid[3:]
        wtres = self.wtres.lower()
        resn = self.resn

        print("Perfroming minimization for mutation {} {} {} ".format(self.wtres, mut_chid_caps, self.resn))
        command = f'mpirun -n {self.thread_num} {self.charmm_dir} -i step3_{wtres}{resn}{mut_chid}_inte.inp > step3_inte_{wtres}{resn}{mut_chid}.out'
        inter_ener_dat_path = "inter_ener_{}{}{}.dat".format(self.wtres, self.resn, mut_chid).lower()
        total_lines = 20
        progress_thread = threading.Thread(target=self.print_progress, args=(inter_ener_dat_path, total_lines))
        progress_thread.start()

        try:
            result = subprocess.run(command, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        except subprocess.CalledProcessError as e:
            print(f"\nError occurred during Step 3: subprocess returned non-zero exit status {e.returncode}.")
            print(f"Python error message: {str(e)}")
            print(f"Please check step3_{wtres}{resn}{mut_chid}.out for more information.")
            self.subprocess_failed = True 
            progress_thread.join()
            sys.exit(e.returncode)

        self.print_progress_exit = True
        progress_thread.join()

        print("\n")

    def print_progress(self, inter_ener_dat_path, total_lines):
        self.print_progress_exit = False
        progress = 0
        while progress < total_lines and not self.print_progress_exit and not self.subprocess_failed:
            time.sleep(1)  # Wait 1 second
            line_count = self.count_lines(inter_ener_dat_path)
            progress = min(line_count, total_lines)
            print(f"\rCalculating progress: {progress}/20", end='', flush=True)

        if self.subprocess_failed:
            sys.exit(2)


    def count_lines(self, inter_ener_dat_path):
        line_count = 0
        while not os.path.exists(inter_ener_dat_path):
            time.sleep(1)  # čaka eno sekundo za refresh
        with open(inter_ener_dat_path, 'r') as inter_ener_dat:
            for line in inter_ener_dat:
                if line.startswith("PRO"):
                    line_count += 1
                    if line_count == 20:
                        return line_count
        return line_count


