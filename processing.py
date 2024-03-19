from flask import Flask, jsonify
from flask_cors import CORS
import subprocess

app = Flask(__name__)
CORS(app)  # Enable CORS for all routes

@app.route('/run-script', methods=['POST'])
def run_script():
    try:
        result = subprocess.run(
            ["python", "C:/Users/acale/OneDrive/Documents/Waterloo BME/4B/BME 462/Interface/sentry.github.io/run_algorithm.py"],
            capture_output=True,
            text=True,
            check=True  # This will raise an exception if the command fails
        )
        return jsonify({"output": result.stdout})
    except subprocess.CalledProcessError as e:
        # Handle the error case
        return jsonify({"error": str(e), "output": e.stdout, "stderr": e.stderr}), 500

if __name__ == "__main__":
    app.run(debug=True)