function evt2ns(evtSamps, peakSamps, varargin)
% EVT2NS Export events to NeuroScope .res / .clu files.
%
%   EVT2NS(evtSamps, peakSamps, varargin)
%
%   SUMMARY:
%       Writes detected events to the NeuroScope format for external
%       visualization, modality-agnostic (ripples, ED, ...). Creates two files
%       tagged by 'fileTag':
%           <basename>.<fileTag>.res.1 - event sample indices (sorted).
%           <basename>.<fileTag>.clu.1 - matching cluster ids.
%       Each event contributes three marks: start (1), end (2), peak (3). Not
%       wired into any pipeline; call it manually to inspect a detection.
%
%   INPUTS:
%       evtSamps    - (Mat)  [N x 2] Event start / end indices (samples).
%       peakSamps   - (Vec)  [N x 1] Event peak indices (samples).
%       varargin    - Parameter/Value pairs:
%           'basepath' - (Char) Target directory. {pwd}
%           'fileTag'  - (Char) Tag between basename and res/clu. {'evt'}
%
%   OUTPUTS:
%       None (writes two files to disk).
%
%   DEPENDENCIES:
%       None.
%
%   HISTORY:
%       Created: 260706 (generalised from ripp2ns; fileTag param; no cd).

%% ========================================================================
%  ARGUMENTS
%  ========================================================================
p = inputParser;
addRequired(p, 'evtSamps', @isnumeric);
addRequired(p, 'peakSamps', @isnumeric);
addParameter(p, 'basepath', pwd, @ischar);
addParameter(p, 'fileTag', 'evt', @ischar);
parse(p, evtSamps, peakSamps, varargin{:});
basepath = p.Results.basepath;
fileTag  = p.Results.fileTag;

[~, basename] = fileparts(basepath);
resfile = fullfile(basepath, [basename, '.', fileTag, '.res.1']);
clufile = fullfile(basepath, [basename, '.', fileTag, '.clu.1']);
nEvents = numel(peakSamps);

%% ========================================================================
%  GENERATE FILES
%  ========================================================================
% Stack start / end / peak samples into one sorted timeline; the cluster id
% (1 = start, 2 = end, 3 = peak) follows the same sort order.
allSamps = [evtSamps(:, 1); evtSamps(:, 2); peakSamps(:)];
[res, sortIdx] = sort(allSamps);
ids = [ones(nEvents, 1); 2 * ones(nEvents, 1); 3 * ones(nEvents, 1)];
clu = ids(sortIdx);

% .res: one sample index per line
fid = fopen(resfile, 'w');
if fid == -1, error('evt2ns:resOpen', 'failed to open %s', resfile); end
fprintf(fid, '%d\n', res);
fclose(fid);

% .clu: leading cluster count, then one id per line
fid = fopen(clufile, 'w');
if fid == -1, error('evt2ns:cluOpen', 'failed to open %s', clufile); end
fprintf(fid, '%d\n', 3);
fprintf(fid, '%d\n', clu);
fclose(fid);

end     % EOF
