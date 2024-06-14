import matlab.engine
import sys
import json

def main():
    try:
        # Start MATLAB engine
        eng = matlab.engine.start_matlab()

        # Add the MATLAB file's path to MATLAB's search path
        eng.addpath("C:/Users/acale/OneDrive/Documents/Waterloo BME/4B/BME 462/Interface/sentry.github.io/")

        print("Current directory:", eng.pwd())
        
        # print(eng.path())

        # Retrieve file path from command line argument
        if len(sys.argv) < 2:
            raise ValueError("File path not provided")

        fpath = sys.argv[1]
        # fpath = "C:/Users/acale/OneDrive/Documents/Waterloo BME/4B/BME 462/CapstoneProjectBase-physician_ui/data/mar13-2024/mar13/"

        # Call the MATLAB function
        result = eng.fullAlgorithm(fpath, nargout=0)
        # print("Algorithm completed successfully")

        # Close MATLAB engine
        eng.quit()

    except Exception as e:
        print("Error occurred:", str(e))
        sys.exit(1)

if __name__ == "__main__":
    main()