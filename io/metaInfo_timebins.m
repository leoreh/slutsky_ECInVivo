function [timebins, timepnt] = metaInfo_timebins(varargin)

% separates the recording length to timebins in relation to specific time
% points determined by the closest value of datInfo.nsamps to the user input
% timepnts.
% 
% INPUT
%   basepath    recording session {pwd}
%   reqPnt      numeric. time of interest from the start of the recording
%               [s]
%   nbins       numeric. number of time bins to create {8}
%   nchans      numeric. number of channels in dat file. if empty will be
%               extracted from session info. used for converting recording
%               length in sample to seconds (same as fs)
%   fs          numeric. sampling frequency of dat file
% 
% OUTPUT
%   timebins    n x 2 mat [s]
%   timepnt     precise point of interest from datInfo.nsamps

p = inputParser;
addOptional(p, 'basepath', pwd);
addOptional(p, 'reqPnt', [], @isnumeric);
addOptional(p, 'nbins', 8, @isnumeric);
addOptional(p, 'fs', [], @isnumeric);
addOptional(p, 'nchans', [], @isnumeric);

parse(p, varargin{:})
basepath        = p.Results.basepath;
reqPnt          = p.Results.reqPnt;
nbins           = p.Results.nbins;
fs              = p.Results.fs;
nchans          = p.Results.nchans;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~, basename] = fileparts(basepath);
sessionfile = fullfile(basepath, [basename, '.session.mat']);
load(sessionfile, 'session')
if isempty(fs)
    fs = session.extracellular.sr;
end
if isempty(nchans)
    nchans = session.extracellular.nChannels;
end
    
% get datInfo file
infofiles = dir('*datInfo*');
fileidx = contains({infofiles.name}, 'EMG');
filename = infofiles(~fileidx).name;
infofile = fullfile(basepath, filename);
load(infofile, 'datInfo')

% recording length
fileinfo = dir([basename, '.dat']);
if isempty(fileinfo)
    fileinfo = dir([basename, '.lfp']);
    recLen = floor(fileinfo.bytes / 2 / nchans / 1250);
else
    recLen = floor(fileinfo.bytes / 2 / nchans / fs);
end
recLen = recLen - 1;    % margin of error

% point of interest relative to 'block' transition in the recording
timepnt = 0;
if ~isempty(reqPnt)
    csec = floor(cumsum(datInfo.nsamps / fs));
    [~, pntIdx] = min(abs(csec - reqPnt));
    timepnt = csec(pntIdx);
end

% timebins
timebins = n2chunks('n', recLen, 'nchunks', nbins, 'pnts', timepnt);

% save session
session.general.timebins = timebins;
session.general.timepnt = timepnt;
save(sessionfile, 'session')

end

% EOF