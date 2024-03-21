document.addEventListener('DOMContentLoaded', function () {
    // Event listener for the 'Process' button
    var runScriptBtn = document.getElementById('run-script-btn');
    if (runScriptBtn) {
        runScriptBtn.addEventListener('click', function() {
            const filePath = 
            "C:\\Users\\acale\\OneDrive\\Documents\\Waterloo BME\\4B\\BME 462\\CapstoneProjectBase-physician_ui\\data\\mar13-2024\\mar13\\";

            // const filePath =
            // "C:\\Users\\acale\\OneDrive\\Documents\\Waterloo BME\\4B\\BME 462\\CapstoneProjectBase-physician_ui\\data\\RBDPatient\\";
            
            fetch('http://127.0.0.1:5000/run-script', { 
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json'
                },
                body: JSON.stringify({filePath: filePath})
            })
            .then(response => {
                if (!response.ok) {
                    throw new Error('Network response was not ok');
                }
                return response.json();
            })
            .then(data => {
                console.log('Script output:', data.output);
                // Notify user
                alert("Processing has been completed successfully. Check the directory for output.");
            })
            .catch(error => {
                console.error('Error:', error);
            });
        });
    }

    // Function to plot EEG from CSV
    function plotEEGFromCSV(eegData) {
        var rows = eegData.split('\n');

        var columnNames = ['Time', 'Delta', 'Theta', 'Alpha', 'Beta', 'Gamma'];
        var xValues = [];
        var yValues = {
            'Delta': [],
            'Theta': [],
            'Alpha': [],
            'Beta': [],
            'Gamma': []
        };
        var yValuesREM = {
            'Delta': [],
            'Theta': [],
            'Alpha': [],
            'Beta': [],
            'Gamma': []
        };

        rows.forEach(function (row) {
            var columns = row.split(',');
            xValues.push(parseFloat(columns[0]));
            for (var i = 1; i <= 5; i++) {
                var yValue = parseFloat(columns[i]);
                yValues[columnNames[i]].push(yValue);

                var remFlag = parseFloat(columns[6]);
                var yValueREM = remFlag !== 0 ? yValue * remFlag : null;
                yValuesREM[columnNames[i]].push(yValueREM);
            }
        });

        var firstXValue = xValues[0];
        xValues = xValues.map(function(value) {
            return value - firstXValue;
        });

        var traceDelta = {
            x: xValues,
            y: yValues[columnNames[1]],
            type: 'scatter',
            name: columnNames[1],
            line: {color: '#4F4698'}
        };
        var traceDeltaREM = {
            x: xValues,
            y: yValuesREM[columnNames[1]],
            type: 'scatter',
            name: columnNames[1],
            line: {color: '#FFB627'}
        };

        var traceTheta = {
            x: xValues,
            y: yValues[columnNames[2]],
            // xaxis: 'x2',
            yaxis: 'y2',
            type: 'scatter',
            name: columnNames[2],
            line: {color: '#4F4698'}
        };
        var traceThetaREM = {
            x: xValues,
            y: yValuesREM[columnNames[2]],
            // xaxis: 'x2',
            yaxis: 'y2',
            type: 'scatter',
            name: columnNames[2],
            line: {color: '#FFB627'}
        };

        var traceAlpha = {
            x: xValues,
            y: yValues[columnNames[3]],
            // xaxis: 'x3',
            yaxis: 'y3',
            type: 'scatter',
            name: columnNames[3],
            line: {color: '#4F4698'}
        };
        var traceAlphaREM = {
            x: xValues,
            y: yValuesREM[columnNames[3]],
            // xaxis: 'x3',
            yaxis: 'y3',
            type: 'scatter',
            name: columnNames[3],
            line: {color: '#FFB627'}
        };

        var traceBeta = {
            x: xValues,
            y: yValues[columnNames[4]],
            // xaxis: 'x4',
            yaxis: 'y4',
            type: 'scatter',
            name: columnNames[4],
            line: {color: '#4F4698'}
        };
        var traceBetaREM = {
            x: xValues,
            y: yValuesREM[columnNames[4]],
            // xaxis: 'x4',
            yaxis: 'y4',
            type: 'scatter',
            name: columnNames[4],
            line: {color: '#FFB627'}
        };

        var traceGamma = {
            x: xValues,
            y: yValues[columnNames[5]],
            // xaxis: 'x5',
            yaxis: 'y5',
            type: 'scatter',
            name: columnNames[5],
            line: {color: '#4F4698'}
        };
        var traceGammaREM = {
            x: xValues,
            y: yValuesREM[columnNames[5]],
            // xaxis: 'x5',
            yaxis: 'y5',
            type: 'scatter',
            name: columnNames[5],
            line: {color: '#FFB627'}
        };

        var eegDataToPlot = [traceGamma, traceGammaREM, traceBeta, traceBetaREM, traceAlpha, traceAlphaREM, traceTheta, traceThetaREM, traceDelta, traceDeltaREM];

        var layout = {
            showlegend: false,
            yaxis: { domain: [0.8, 1], title: 'Delta', fixedrange: true },
            yaxis2: { domain: [0.6, 0.79], title: 'Theta', fixedrange: true },
            yaxis3: { domain: [0.4, 0.59], title: 'Alpha', fixedrange: true },
            yaxis4: { domain: [0.2, 0.39], title: 'Beta', fixedrange: true },
            yaxis5: { domain: [0, 0.19], title: 'Gamma', fixedrange: true },
            xaxis: { title: 'Time (seconds)', fixedrange: false, anchor: 'y5', domain: [0, 1] },
            dragmode: 'pan',
            zoom: {
                type: 'x',
                wheel: 'x',
                dragmode: 'zoom'
            }
        };

        Plotly.newPlot('eegPlot', eegDataToPlot, layout);

    }

    // Function to plot IMU from CSV
    function plotIMUFromCSV(imuData) {
        var colours = ['#4F4698', '#FFD166', '#CC92C2', '#FFB627'];
        var traces = [];

        for (const [index, key] of Object.keys(imuData).entries()) {
            var data = imuData[key];
            var xValues = data.map(row => parseFloat(row[0]));
            var yValues = data.map(row => parseFloat(row[1]));

            var firstXValue = xValues[0];
            xValues = xValues.map(function(value) {
                return value - firstXValue;
            });

            var trace = {
                x: xValues,
                y: yValues,
                mode: 'lines',
                name: key,
                line: {color: colours[index]}
            };

            traces.push(trace);
        }

        var layout = {
            xaxis: {title: 'Time (seconds)'},
            yaxis: {title: 'Magnitude'},
            legend: {
                x: 0,
                y: 1,
                xanchor: 'bottom',
                yanchor: 'bottom'
            }
        };

        Plotly.newPlot('imuPlot', traces, layout);
    }

    // Event listener for file input change
    document.getElementById('inputFiles').addEventListener('change', function (event) {
        const allFiles = event.target.files;
        if (!allFiles || allFiles.length === 0) return;

        const eegFile = Array.from(allFiles).find(file => file.name === "eeg.csv");

        if (!eegFile) {
            console.error("EEG file 'eeg.csv' not found.");
            return;
        }

        const imuFiles = Array.from(allFiles).filter(file => {
            return file.name === "leftWrist.csv" || file.name === "rightWrist.csv" || file.name === "leftAnkle.csv" || file.name === "rightAnkle.csv";
        });

        const eegReader = new FileReader();

        eegReader.onload = function (e) {
            const eegData = e.target.result;
            plotEEGFromCSV(eegData);
        };

        eegReader.readAsText(eegFile);

        const imuData = {};
        let imuFilesCount = 0;

        imuFiles.forEach(file => {
            const imuReader = new FileReader();
            imuReader.onload = function(e) {
                const csvData = e.target.result;
                const csvLines = csvData.split('\n');
                const dataArray = csvLines.map(line => line.split(','));
                const fileName = file.name.split('.')[0];
                imuData[fileName] = dataArray;
                imuFilesCount++;

                if (imuFilesCount === imuFiles.length) {
                    plotIMUFromCSV(imuData);
                    console.log(imuData);
                }
            };
            imuReader.readAsText(file);
        });
    });
});