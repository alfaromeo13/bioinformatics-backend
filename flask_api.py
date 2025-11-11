import os
import uuid
import zipfile
import subprocess
from flask_cors import CORS
from collections import deque
from flask import Flask, Response, request, jsonify, send_file, abort

app = Flask(__name__) # Flask constructor.

'''
https://stackoverflow.com/questions/25594893/how-to-enable-cors-in-flask
This will only allow CORS requests from http://localhost:4200
'''
CORS(app, resources={r"/*": {"origins": "http://localhost:4200"}})
app.config['CORS_HEADERS'] = 'Content-Type'

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_FOLDER = os.path.join(BASE_DIR, 'outputs')
os.makedirs(OUTPUT_FOLDER, exist_ok=True)

@app.route('/run-script', methods=['POST'])
def run_script():
    try:
        pdb_file = request.files.get('pdb_file')
        if not pdb_file:
            return jsonify({"error": "PDB file is required"}), 400

        protein_chains = request.form.get('protein_chains', '').strip()
        partner_chains = request.form.get('partner_chains', '').strip()
        mutations = request.form.get('mutations', '').strip()
        detect_interface = request.form.get('detect_interface', 'false') == 'true'

        job_id = str(uuid.uuid4())[:8]
        job_prefix = f"job_{job_id}"
        pdb_name = f"{job_prefix}_{pdb_file.filename}"
        pdb_path = os.path.join(BASE_DIR, pdb_name)
        pdb_file.save(pdb_path)
        print(f"[IIME] Saved {pdb_path}")

        container_name = f"iime_{job_id}"

        docker_cmd = [
            "docker", "run", "--name", container_name,
            "-v", f"{BASE_DIR}:/IIME",
            "-w", "/IIME",
            "iime_env_full",
            "python3", "IIME.py",
            "--charmm-dir", "/usr/local/charmm_program/charmm/bin/charmm",
            "--pdb", pdb_name,
            "--protein-chains", protein_chains,
            "--partner-chains", partner_chains,
            "-l", f"run_{job_id}.log",
            "--threads", "4"
        ]

        if mutations and detect_interface:
            return jsonify({
                "error": "Mutations and Detect Interface are mutually exclusive. Please select only one."
            }), 400

        if mutations:
            docker_cmd += ["--mutations", mutations]
        elif detect_interface:
            docker_cmd += ["--cutoff", "5.0"]

        subprocess.Popen(docker_cmd, cwd=BASE_DIR)
        return jsonify({"status": "processing", "job_id": job_id})

    except Exception as e:
        return jsonify({"error": str(e)}), 500

@app.route('/check-result/<job_id>', methods=['GET'])
def check_result(job_id):
    """Check results specific to a job"""
    zips = [f for f in os.listdir(BASE_DIR) if f.endswith('_results.zip') and job_id in f]
    if not zips:
        return jsonify({"status": "pending"}), 202

    zip_path = os.path.join(BASE_DIR, zips[0])
    extract_dir = os.path.join(OUTPUT_FOLDER, job_id)

    if not os.path.exists(extract_dir) or not os.listdir(extract_dir):
        os.makedirs(extract_dir, exist_ok=True)
        with zipfile.ZipFile(zip_path, 'r') as zip_ref:
            zip_ref.extractall(extract_dir)
        print(f"[IIME] Extracted to {extract_dir}")
    else:
        print(f"[IIME] {job_id} already extracted.")

    files = []
    for root, _, filenames in os.walk(extract_dir):
        for fname in filenames:
            files.append(os.path.relpath(os.path.join(root, fname), extract_dir))

    return jsonify({"status": "completed", "files": files})

@app.route('/get-file/<job_id>/<path:filename>', methods=['GET'])
def get_file(job_id, filename):
    """Return extracted result file (text or binary)"""
    extract_dir = os.path.join(OUTPUT_FOLDER, job_id)
    file_path = os.path.join(extract_dir, filename)

    if os.path.exists(file_path):
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                content = f.read()
            return Response(content, mimetype='text/plain')
        except UnicodeDecodeError:
            return send_file(file_path, as_attachment=True)

    abort(404, description=f"File not found: {filename}")

@app.route('/get-log/<job_id>', methods=['GET'])
def get_log(job_id):
    container_name = f"iime_{job_id}"
    try:
        result = subprocess.run(
            ["docker", "logs", "--tail", "10", container_name],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            timeout=5
        )
        return jsonify({"log": result.stdout}), 200
    except subprocess.CalledProcessError:
        return jsonify({"log": "[Container finished or not found]"}), 200


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