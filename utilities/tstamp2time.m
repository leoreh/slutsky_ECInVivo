function [dt, tstamp] = tstamp2time(varargin)

% converts timestamp to datetime based on the input date time string.
% designed to work with basename (mousename_yyMMdd_HHmmss). for example, if
% dtstr = 210228_190000 and tstamp = 7200 then dt = 28 feb 21 21:00:00. can
% also work the otherway around; if tstamp is empty and dtstr is different
% than basename, will return the second within the recording that
% corresponds to dtstr. if sampling frequency is specified, than assumes
% tstamp is in samples rather than seconds and will also return dt in
% samples

% INPUT:
%   basepath            string. path to recording folder {pwd}.
%   dtstr               string of date time or datetime. can also be basename in the
%                       format "mousename_date_time"
%   dtFormat            string. format of dtstr. if empty will be extracted
%                       from length of dtstr. best is 'yyMMdd_HHmmss'
%   tstr                string of time in the same time format as dtstr. 
%   tstamp              numeric. seconds / samples to add (positive) or remove
%                       (negative) from dtstr [seconds]
%   fs                  sampling frequency. if specified will assume
%                       tstamp is in samples
% 
% DEPENDENCIES
%   guessDateTime
%
% TO DO LIST
%   # allow tstamp to be in samples by inputing fs (done)
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
addOptional(p, 'tstamp', [], @isnumeric);
addOptional(p, 'fs', 1, @isnumeric);

parse(p, varargin{:})
dtstr               = p.Results.dtstr;
dtFormat            = p.Results.dtFormat;
tstr                = p.Results.tstr;
tstamp              = p.Results.tstamp;
fs                  = p.Results.fs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% conversion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isdatetime(dtstr)
    [dtBase, ~] = guessDateTime(dtstr);
else
    dtBase = dtstr;
end

% ALT 1: dtstr and tstamp specified; will add / remove tstamp from dtstr
% and return output as date time
if ~isempty(tstamp)
    dt = dtBase + seconds(tstamp / fs);
else
    dt = [];
end

% ALT 2: dtstr and tstr specified; will find tstamp of tstr relative to
% dtstr in samples / seconds
if ~isempty(tstr)
    if ~isdatetime(tstr)
        [dtTime, ~] = guessDateTime(tstr);
    else
        dtTime = tstr;
    end
    tstamp(1) = seconds(dtTime - dtBase) * fs;
%     tstamp(2) = seconds(dtTime + hours(12) - dtBase) * fs;
end

end

