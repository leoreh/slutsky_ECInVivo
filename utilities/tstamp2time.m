function [dt, tstamps] = tstamps2time(varargin)

% converts timestamp to datetime based on the input date time string.
% designed to work with basename (mousename_yyMMdd_HHmmss). for example, if
% dtstr = 210228_190000 and tstamps = 7200 then dt = 28 feb 21 21:00:00. can
% also work the otherway around; if tstamps is empty and dtstr is different
% than basename, will return the second within the recording that
% corresponds to dtstr. if sampling frequency is specified, than assumes
% tstamps is in samples rather than seconds and will also return dt in
% samples

% INPUT:
%   basepath            string. path to recording folder {pwd}.
%   dtstr               string of date time or datetime. can also be basename in the
%                       format "mousename_date_time"
%   dtFormat            string. format of dtstr. if empty will be extracted
%                       from length of dtstr. best is 'yyMMdd_HHmmss'
%   tstr                string of time in the same time format as dtstr. 
%   tstamps              numeric. seconds / samples to add (positive) or remove
%                       (negative) from dtstr [seconds]
%   fs                  sampling frequency. if specified will assume
%                       tstamps is in samples
% 
% DEPENDENCIES
%   guessDateTime
%
% TO DO LIST
%   # allow tstamps to be in samples by inputing fs (done)
% 
% 02 apr 21 LH  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;

p = inputParser;
addOptional(p, 'dtstr', []);
addOptional(p, 'dtFormat', []);
addOptional(p, 'tstr', []);
addOptional(p, 'tstamps', [], @isnumeric);
addOptional(p, 'fs', 1, @isnumeric);

parse(p, varargin{:})
dtstr               = p.Results.dtstr;
dtFormat            = p.Results.dtFormat;
tstr                = p.Results.tstr;
tstamps              = p.Results.tstamps;
fs                  = p.Results.fs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% conversion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isdatetime(dtstr)
    [dtBase, ~] = guessDateTime(dtstr);
else
    dtBase = dtstr;
end

% ALT 1: dtstr and tstamps specified; will add / remove tstamps from dtstr
% and return output as date time
if ~isempty(tstamps)
    for it = 1 : length(tstamps)
        dt = dtBase + seconds(tstamps / fs);
    end
else
    dt = [];
end

% ALT 2: dtstr and tstr specified; will find tstamps of tstr relative to
% dtstr in samples / seconds
if ~isempty(tstr)
    if ~isdatetime(tstr)
        [dtTime, ~] = guessDateTime(tstr);
    else
        dtTime = tstr;
    end
    tstamps(1) = seconds(dtTime - dtBase) * fs;
end

end

