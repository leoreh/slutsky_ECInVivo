function ripp2ns(rippSamps, peakSamps, varargin)
% RIPP2NS Generates NeuroScope-compatible .res and .clu files.
%
%   RIPP2NS(rippSamps, peakSamps, varargin)
%
%   SUMMARY:
%       Exports ripple events to the NeuroScope format for external visualization.
%       Creates two files:
%           - .res.X: Event timestamps (in samples).
%           - .clu.X: Event Cluster IDs.
%
%       Mapping:
%           Cluster 1: Ripple Start
%           Cluster 2: Ripple End
%           Cluster 3: Ripple Peak
%
%   INPUTS:
%       rippSamps   - (Mat) [N x 2] Start/End indices (Sample count).
%       peakSamps   - (Vec) [N x 1] Peak indices (Sample count).
%       varargin    - Parameter/Value pairs:
%           'basepath' - (Char) Target directory. (Default: pwd).
%
%   OUTPUTS:
%       None. Saves files to disk.
%
%   DEPENDENCIES:
%       None.
%
%   HISTORY:
%       Updated: 23 Jan 2026

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



end     % EOF
