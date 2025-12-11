function mcu_dashboard(frTbl, tAxis, varargin)

% MCU_DASHBOARD Creates an interactive HTML GUI for data exploration.
%   Generates a standalone HTML file with interactive Plotly graphs.
%   Users can filter by Name and Unit Type, grouped by 'Group'.
%   Includes CSV download and refined visuals (Mean on top, Black).
%
% INPUTS:
%   frTbl     (table) - Table containing 'FRt', 'Name', 'UnitType', and optional 'Group'.
%   tAxis     (vec)   - Time axis vector corresponding to columns of FRt.
%   filename  (char)  - Optional. Name of the output HTML file.
%                       Default: 'Data_Dashboard.html'
%
% OUTPUTS:
%   None. Saves a file to disk.
%
% See also: MCU_FRTBL

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 'frTbl', @istable);
addRequired(p, 'tAxis', @isnumeric);
addOptional(p, 'filename', 'Data_Dashboard.html', @ischar);
parse(p, frTbl, tAxis, varargin{:});

filename = p.Results.filename;

%% ========================================================================
%  PREPARE DATA
%  ========================================================================

% Remove rows with NO Firing Rate data (all NaNs)
validIdx = any(~isnan(frTbl.FRt), 2);
cleanTbl = frTbl(validIdx, :);

% Ensure Group column exists
if ~ismember('Group', cleanTbl.Properties.VariableNames)
    cleanTbl.Group = repmat({'All'}, height(cleanTbl), 1);
end

% 1. The Time Axis
jsonTime = jsonencode(tAxis);

% 2. The Metadata (Array of objects: [{m:'lh132', t:'RS', g:'WT'}, ...])
% Convert categorical to string for JSON safety
if iscategorical(cleanTbl.Name), cleanTbl.Name = cellstr(cleanTbl.Name); end
if iscategorical(cleanTbl.UnitType), cleanTbl.UnitType = cellstr(cleanTbl.UnitType); end
if iscategorical(cleanTbl.Group), cleanTbl.Group = cellstr(cleanTbl.Group); end

% Usage of 'Name' key for mouse identifier
metaStruct = table2struct(cleanTbl(:, {'Name', 'UnitType', 'Group'}));

% Rename fields to be short (minimize file size)
[metaStruct.m] = metaStruct.Name;
[metaStruct.u] = metaStruct.UnitType;
[metaStruct.g] = metaStruct.Group;
metaStruct = rmfield(metaStruct, {'Name', 'UnitType', 'Group'});
jsonMeta = jsonencode(metaStruct);

% 3. The Firing Rate Matrix (Array of arrays)
% Round to 4 decimals to reduce file size significantly
frMatrix = round(cleanTbl.FRt, 4);
jsonData = jsonencode(frMatrix);

%% ========================================================================
%  WRITE HTML FRAMEWORK
%  ========================================================================

fid = fopen(filename, 'w', 'n', 'utf-8');

fprintf(fid, '<html><head><meta charset="utf-8"/><title>Perturbation Dashboard</title>');
fprintf(fid, '<script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>');
fprintf(fid, '<style>');
fprintf(fid, '  body { font-family: "Segoe UI", sans-serif; display: flex; height: 100vh; margin: 0; }');
fprintf(fid, '  #sidebar { width: 250px; background: #f4f4f4; padding: 20px; border-right: 1px solid #ccc; overflow-y: auto; flex-shrink: 0; }');
fprintf(fid, '  #main { flex-grow: 1; padding: 20px; display: flex; flex-direction: column; overflow: hidden; }');
fprintf(fid, '  #plotDiv { flex-grow: 1; min-height: 400px; }');
fprintf(fid, '  .control-group { margin-bottom: 20px; }');
fprintf(fid, '  label { display: block; margin: 5px 0; cursor: pointer; }');
fprintf(fid, '  h3 { margin-top: 0; font-size: 16px; color: #333; border-bottom: 2px solid #ccc; padding-bottom: 5px; }');
fprintf(fid, '  h4 { margin: 10px 0 5px; font-size: 14px; color: #555; text-transform: uppercase; font-weight: bold; }');
fprintf(fid, '  button { padding: 10px 20px; background: #0072BD; color: white; border: none; cursor: pointer; width: 100%%; margin-bottom: 10px; }');
fprintf(fid, '  button:hover { background: #005bd3; }');
fprintf(fid, '  .btn-green { background: #28a745; }');
fprintf(fid, '  .btn-green:hover { background: #218838; }');
fprintf(fid, '</style></head><body>');

%% ========================================================================
%  SIDEBAR CONTROLS
%  ========================================================================

fprintf(fid, '<div id="sidebar">');
fprintf(fid, '  <h2>Data Controls</h2>');

% Unit Type Checkboxes
fprintf(fid, '  <div class="control-group"><h3>Unit Type</h3>');
uTypes = unique(cleanTbl.UnitType);
for i = 1:length(uTypes)
    uVal = char(uTypes(i));
    fprintf(fid, '<label><input type="checkbox" class="u-check" value="%s" checked> %s</label>', uVal, uVal);
end
fprintf(fid, '  </div>');

% Mouse Checkboxes - GROUPED
fprintf(fid, '  <div class="control-group"><h3>Mice</h3>');
groups = unique(cleanTbl.Group);

% If we have Groups, sort/display by Group
for iG = 1:length(groups)
    grpName = char(groups(iG));

    % Find mice in this group
    grpMice = unique(cleanTbl.Name(strcmp(cleanTbl.Group, grpName)));

    if ~isempty(grpMice)
        fprintf(fid, '  <h4>%s</h4>', grpName);
        for iM = 1:length(grpMice)
            mVal = char(grpMice(iM));
            fprintf(fid, '<label><input type="checkbox" class="m-check" value="%s" checked> %s</label>', mVal, mVal);
        end
    end
end
fprintf(fid, '  </div>');

% Buttons
fprintf(fid, '  <button onclick="updatePlot()">Update Graph</button>');
fprintf(fid, '  <button class="btn-green" onclick="downloadCSV()">Download CSV</button>');
fprintf(fid, '  <br><br><small id="status">Ready</small>');
fprintf(fid, '</div>');

%% ========================================================================
%  MAIN PLOT AREA
%  ========================================================================

fprintf(fid, '<div id="main">');
fprintf(fid, '  <h1>Interactive Firing Rate Explorer</h1>');
fprintf(fid, '  <div id="plotDiv"></div>');
fprintf(fid, '</div>');

%% ========================================================================
%  JAVASCRIPT LOGIC
%  ========================================================================

fprintf(fid, '<script>');

% Embed the Data
fprintf(fid, 'const time = %s;', jsonTime);
fprintf(fid, 'const meta = %s;', jsonMeta);
fprintf(fid, 'const frData = %s;', jsonData);

% --- Helper: Get Filtered Indices ---
fprintf(fid, 'function getFilteredIndices() {');
fprintf(fid, '  var selUnits = Array.from(document.querySelectorAll(".u-check:checked")).map(cb => cb.value);');
fprintf(fid, '  var selMice = Array.from(document.querySelectorAll(".m-check:checked")).map(cb => cb.value);');
fprintf(fid, '  var indices = [];');
fprintf(fid, '  for (var i = 0; i < meta.length; i++) {');
fprintf(fid, '    if (selUnits.includes(meta[i].u) && selMice.includes(meta[i].m)) {');
fprintf(fid, '      indices.push(i);');
fprintf(fid, '    }');
fprintf(fid, '  }');
fprintf(fid, '  return indices;');
fprintf(fid, '}');

% --- Main Plotting Function ---
fprintf(fid, 'function updatePlot() {');
fprintf(fid, '  var status = document.getElementById("status");');
fprintf(fid, '  status.innerText = "Processing...";');
fprintf(fid, '  var indices = getFilteredIndices();');

fprintf(fid, '  if (indices.length === 0) { Plotly.purge("plotDiv"); status.innerText = "No data selected"; return; }');

% 3. Aggregate Data for Plotting
fprintf(fid, '  var traces = [];');

% Calculate Mean
fprintf(fid, '  var meanTrace = { x: time, y: [], mode: "lines", line: {color: "black", width: 4}, name: "MEAN (" + indices.length + " units)" };');
fprintf(fid, '  for (var t = 0; t < time.length; t++) {');
fprintf(fid, '    var sum = 0; var count = 0;');
fprintf(fid, '    for (var k = 0; k < indices.length; k++) {');
fprintf(fid, '       var val = frData[indices[k]][t];');
fprintf(fid, '       if (val !== null && !isNaN(val)) { sum += val; count++; }');
fprintf(fid, '    }');
fprintf(fid, '    meanTrace.y.push(count > 0 ? sum/count : null);');
fprintf(fid, '  }');

% Add Sampled Individual Traces FIRST (so they are behind the mean)
fprintf(fid, '  var maxTraces = 150;');
fprintf(fid, '  var step = Math.ceil(indices.length / maxTraces);');
fprintf(fid, '  for (var k = 0; k < indices.length; k += step) {');
fprintf(fid, '     var idx = indices[k];');
fprintf(fid, '     traces.push({');
fprintf(fid, '       x: time, y: frData[idx],');
fprintf(fid, '       mode: "lines", line: {color: "rgba(150,150,150,0.3)", width: 1},');
fprintf(fid, '       hoverinfo: "skip", showlegend: false');
fprintf(fid, '     });');
fprintf(fid, '  }');

% Push Mean Trace LAST (so it renders on top)
fprintf(fid, '  traces.push(meanTrace);');

% 4. Render
fprintf(fid, '  var layout = { title: "", xaxis: {title: "Time (h)", tickmode: "linear", tick0: 0, dtick: 24}, yaxis: {title: "Firing Rate (Hz)"} };');
fprintf(fid, '  Plotly.newPlot("plotDiv", traces, layout);');
fprintf(fid, '  status.innerText = "Done. Showing " + indices.length + " units.";');
fprintf(fid, '}');

% --- CSV Download Function ---
fprintf(fid, 'function downloadCSV() {');
fprintf(fid, '  var status = document.getElementById("status");');
fprintf(fid, '  status.innerText = "Generating CSV...";');
fprintf(fid, '  var indices = getFilteredIndices();');
fprintf(fid, '  if (indices.length === 0) { alert("No data to export"); return; }');
fprintf(fid, '  ');
fprintf(fid, '  var csvContent = "data:text/csv;charset=utf-8,";');
fprintf(fid, '  ');
% Header Row
fprintf(fid, '  var header = ["Time"];');
fprintf(fid, '  indices.forEach(idx => { header.push(meta[idx].m + "_" + meta[idx].u + "_" + idx); });');
fprintf(fid, '  csvContent += header.join(",") + "\\n";');
fprintf(fid, '  ');
% Data Rows
fprintf(fid, '  for (var t = 0; t < time.length; t++) {');
fprintf(fid, '      var row = [time[t]];');
fprintf(fid, '      indices.forEach(idx => { row.push(frData[idx][t]); });');
fprintf(fid, '      csvContent += row.join(",") + "\\n";');
fprintf(fid, '  }');
fprintf(fid, '  ');
% Download Trigger
fprintf(fid, '  var encodedUri = encodeURI(csvContent);');
fprintf(fid, '  var link = document.createElement("a");');
fprintf(fid, '  link.setAttribute("href", encodedUri);');
fprintf(fid, '  link.setAttribute("download", "mcu_data_export.csv");');
fprintf(fid, '  document.body.appendChild(link);');
fprintf(fid, '  link.click();');
fprintf(fid, '  document.body.removeChild(link);');
fprintf(fid, '  status.innerText = "Export Complete.";');
fprintf(fid, '}');


% Initialize
fprintf(fid, 'updatePlot();');
fprintf(fid, '</script></body></html>');

fclose(fid);
fprintf('Dashboard saved: <a href="matlab:web(''%s'')">%s</a>\n', filename, filename);

end