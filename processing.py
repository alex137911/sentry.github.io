from flask import Flask, request, jsonify
from flask_cors import CORS
import subprocess

app = Flask(__name__)
CORS(app)

@app.route('/run-script', methods=['POST'])
def run_script():
    data = request.get_json()
    print("Received data:", data)  # This line logs the received data
    fpath = data.get('filePath')

    if not fpath:
        return jsonify({"error": "File path not provided"}), 400

    try:
        # ABSOLUTE PATH (DOES NOT CHANGE)
        result = subprocess.run(
            ["C:/ProgramData/anaconda3/python.exe", # Path to conda python (3.11.5)
             "C:/Users/acale/OneDrive/Documents/Waterloo BME/4B/BME 462/Interface/sentry.github.io/run_algorithm.py", 
             fpath],
            capture_output=True,
            text=True,
            check=True
        )
        return jsonify({"output": result.stdout})
    
    except subprocess.CalledProcessError as e:
        return jsonify({"error": str(e), "output": e.stdout, "stderr": e.stderr}), 500

if __name__ == "__main__":
    app.run(debug=True)