import os
import uuid
import shutil
import zipfile
import subprocess
from collections import deque
from flask import Flask, Response, request, jsonify, send_file, abort

app = Flask(__name__)

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


@app.route('/run-script', methods=['POST'])
def run_script():

    # Cleanup old outputs
    for folder in os.listdir(OUTPUT_FOLDER):
        shutil.rmtree(os.path.join(OUTPUT_FOLDER, folder), ignore_errors=True)

    try:
        pdb_file = request.files.get('pdb_file')
        if not pdb_file:
            return jsonify({"error": "PDB file is required"}), 400

        protein_chains = request.form.get('protein_chains', '').strip()
        partner_chains = request.form.get('partner_chains', '').strip()
        mutations = request.form.get('mutations', '').strip()

        job_id = str(uuid.uuid4())[:8]
        job_prefix = f"job_{job_id}"

        pdb_name = f"{job_prefix}_{pdb_file.filename}"
        pdb_path = os.path.join(BASE_DIR, pdb_name)
        pdb_file.save(pdb_path)

        log_path = os.path.join(LOG_FOLDER, f"run_{job_id}.log")

        cmd = [
            "python3", "IIME.py",
            "--pdb", pdb_name,
            "--protein-chains", protein_chains,
            "--partner-chains", partner_chains,
            "-l", log_path,
            "--threads", "4"
        ]

        if mutations:
            cmd += ["--mutations", mutations]

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
        job["log_file"].close()
        RUNNING_JOBS.pop(job_id, None)

        try:
            os.remove(log_path)
        except:
            pass

        log_content += "\n[Job finished, log deleted]"

    return jsonify({"log": log_content}), 200


@app.route('/check-result/<job_id>', methods=['GET'])
def check_result(job_id):
    zips = [f for f in os.listdir(BASE_DIR) if f.endswith('_results.zip') and job_id in f]
    if not zips:
        return jsonify({"status": "pending"}), 202

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

    return jsonify({"status": "completed", "files": files})


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
    app.run(host='0.0.0.0', port=3001)