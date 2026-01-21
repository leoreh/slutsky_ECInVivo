function ripp2ns(rippSamps, peakSamps, varargin)
% RIPP2NS Generates NeuroScope .res and .clu files for ripple events.
%
% SUMMARY:
%   Creates .res and .clu files compatible with NeuroScope visualization.
%   Events are stored as: Start (1), End (2), Peak (3).
%
% INPUTS:
%   rippSamps   - (N x 2) Matrix of ripple start and end samples.
%   peakSamps   - (N x 1) Vector of ripple peak samples.
%   varargin    - (param/value) Optional parameters:
%                 'basepath'  : (char) Session path (default: pwd).
%
% OUTPUTS:
%   None. Writes .res and .clu files to disk.
%
% DEPENDENCIES: None

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 'rippSamps', @isnumeric);
addRequired(p, 'peakSamps', @isnumeric);
addParameter(p, 'basepath', pwd, @ischar);

parse(p, rippSamps, peakSamps, varargin{:});
basepath = p.Results.basepath;


%% ========================================================================
%  PREPARATION
%  ========================================================================

cd(basepath);
[~, basename] = fileparts(basepath);
resfile = fullfile(basepath, [basename, '.ripp.res.1']);
clufile = fullfile(basepath, [basename, '.ripp.clu.1']);
nEvents = length(peakSamps);

%% ========================================================================
%  GENERATE FILES
%  ========================================================================

fprintf('Generating NeuroScope files for %d events...\n', nEvents);

% Combine all samples
allSamps = [rippSamps(:, 1); rippSamps(:, 2); peakSamps];

% Sort events by sample index
[res, sort_idx] = sort(allSamps);

% Create CLU ID vector
% 1 = Start, 2 = End, 3 = Peak
ids = [ones(nEvents, 1); ones(nEvents, 1) * 2; ones(nEvents, 1) * 3];
clu = ids(sort_idx);

% Write .res file
fid = fopen(resfile, 'w');
if fid == -1, error('Failed to open .res file.'); end
fprintf(fid, '%d\n', res);
fclose(fid);

% Write .clu file
fid = fopen(clufile, 'w');
if fid == -1, error('Failed to open .clu file.'); end
fprintf(fid, '%d\n', 3); % Number of clusters/classes
fprintf(fid, '%d\n', clu);
fclose(fid);

fprintf('Saved: %s\n', resfile);

end     % EOF
