import os
import sys
import uuid
import shutil
import zipfile
import logging
import subprocess
from flask_cors import CORS
from collections import deque
from flask import Flask, Response, request, jsonify, send_file, abort

# Configure logging to stdout for PM2
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s',
    handlers=[logging.StreamHandler(sys.stdout)]
)
logger = logging.getLogger(__name__)

app = Flask(__name__)

CORS(
    app,
    supports_credentials=True,
    resources={r"/*": {"origins": "*"}}
)

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_FOLDER = os.path.join(BASE_DIR, 'outputs')
LOG_FOLDER = os.path.join(BASE_DIR, 'logs')

os.makedirs(OUTPUT_FOLDER, exist_ok=True)
os.makedirs(LOG_FOLDER, exist_ok=True)

# Track running jobs
RUNNING_JOBS = {}

# -------------------------
# Utility: read last N lines
# -------------------------
def tail_file(path, lines=10):
    try:
        with open(path, "r", encoding="utf-8", errors="ignore") as f:
            return "".join(deque(f, maxlen=lines))
    except FileNotFoundError:
        return "[Log file not found]"



@app.before_request
def handle_options():
    if request.method == "OPTIONS":
        return "", 200


@app.route('/run-script', methods=['POST'])
def run_script():
    logger.info("Received /run-script request")

    # Cleanup old outputs
    for folder in os.listdir(OUTPUT_FOLDER):
        shutil.rmtree(os.path.join(OUTPUT_FOLDER, folder), ignore_errors=True)

    try:
        pdb_file = request.files.get('pdb_file')
        if not pdb_file:
            logger.error("PDB file is required but not provided")
            return jsonify({"error": "PDB file is required"}), 400

        protein_chains = request.form.get('protein_chains', '').strip()
        partner_chains = request.form.get('partner_chains', '').strip()
        mutations = request.form.get('mutations', '').strip()

        job_id = str(uuid.uuid4())[:8]
        job_prefix = f"job_{job_id}"

        logger.info(f"Starting job {job_id}: pdb={pdb_file.filename}, protein_chains={protein_chains}, partner_chains={partner_chains}, mutations={mutations}")

        pdb_name = f"{job_prefix}_{pdb_file.filename}"
        pdb_path = os.path.join(BASE_DIR, pdb_name)
        pdb_file.save(pdb_path)

        log_path = os.path.join(LOG_FOLDER, f"run_{job_id}.log")

        cmd = [
            "python3", "IIME.py",
            "--charmm-dir", os.path.join(BASE_DIR, "charmm_program", "charmm", "bin", "charmm"),
            "--pdb", pdb_name,
            "--protein-chains", protein_chains,
            "--partner-chains", partner_chains,
            "-l", log_path,
            "--threads", "4"
        ]

        if mutations:
            cmd += ["--mutations", mutations]

        logger.info(f"Running command: {' '.join(cmd)}")

        log_file = open(log_path, "w")

        proc = subprocess.Popen(
            cmd,
            cwd=BASE_DIR,
            stdout=log_file,
            stderr=subprocess.STDOUT
        )

        RUNNING_JOBS[job_id] = {
            "process": proc,
            "log_path": log_path,
            "log_file": log_file
        }

        return jsonify({"status": "processing", "job_id": job_id})

    except Exception as e:
        logger.exception(f"Error in /run-script: {e}")
        return jsonify({"error": str(e)}), 500


@app.route('/get-log/<job_id>', methods=['GET'])
def get_log(job_id):
    job = RUNNING_JOBS.get(job_id)

    if not job:
        return jsonify({"log": "[No such job]"}), 200

    proc = job["process"]
    log_path = job["log_path"]

    log_content = tail_file(log_path, lines=10)

    # If process finished → cleanup
    if proc.poll() is not None:
        exit_code = proc.returncode
        job["log_file"].close()
        RUNNING_JOBS.pop(job_id, None)

        if exit_code != 0:
            logger.error(f"Job {job_id} failed with exit code {exit_code}")
            logger.error(f"Last log output:\n{log_content}")
        else:
            logger.info(f"Job {job_id} completed successfully")

        try:
            os.remove(log_path)
        except:
            pass

        log_content += f"\n[Job finished with exit code {exit_code}, log deleted]"

    return jsonify({"log": log_content}), 200


@app.route('/check-result/<job_id>', methods=['GET'])
def check_result(job_id):
    try:
        zips = [f for f in os.listdir(BASE_DIR) if f.endswith('_results.zip') and job_id in f]
        if not zips:
            return jsonify({"status": "pending"}), 202

        logger.info(f"Found results for job {job_id}: {zips[0]}")

        zip_path = os.path.join(BASE_DIR, zips[0])
        extract_dir = os.path.join(OUTPUT_FOLDER, job_id)

        if not os.path.exists(extract_dir):
            os.makedirs(extract_dir, exist_ok=True)
            with zipfile.ZipFile(zip_path, 'r') as zip_ref:
                zip_ref.extractall(extract_dir)

        files = []
        for root, _, filenames in os.walk(extract_dir):
            for fname in filenames:
                files.append(os.path.relpath(os.path.join(root, fname), extract_dir))

        logger.info(f"Job {job_id} completed with {len(files)} result files")
        return jsonify({"status": "completed", "files": files})
    except Exception as e:
        logger.exception(f"Error checking result for job {job_id}: {e}")
        return jsonify({"error": str(e)}), 500


@app.route('/get-file/<job_id>/<path:filename>', methods=['GET'])
def get_file(job_id, filename):
    extract_dir = os.path.join(OUTPUT_FOLDER, job_id)
    file_path = os.path.join(extract_dir, filename)

    if not os.path.exists(file_path):
        abort(404)

    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            return Response(f.read(), mimetype='text/plain')
    except UnicodeDecodeError:
        return send_file(file_path, as_attachment=True)


# Following code is run only when the script is executed, not when it’s imported as a module!
if __name__ == '__main__':
    '''
    The line below starts the Flask development server
    Basically, it runs a web server (built into Flask) on your machine.
    Listens for HTTP requests (like from a browser or frontend).
    Routes requests to the appropriate functions (based on @app.route() decorators).
    Using 0.0.0.0 exposes your app to the entire network or internet (Any computer that can reach servers IP)
    Refer to the following link for more: https://stackoverflow.com/questions/7023052/configure-flask-dev-server-to-be-visible-across-the-network
    '''
    app.run(host='0.0.0.0', port=5000)