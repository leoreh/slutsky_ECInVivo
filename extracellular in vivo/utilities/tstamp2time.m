function dt = tstamp2time(varargin)

% converts timestamp to datetime based on the input date time string.
% designed to work with basename (mousename_yyMMdd_HHmmss). for example, if
% dtstr = 210228_190000 and tstamp = 7200 than dt = 28 feb 21 21:00:00. can
% also work the otherway around; if tstamp is empty and dtstr is different
% than basename, will return the second within the recording the
% corresponds to dtstr. if sampling frequency is specified, than assumes
% tstamp is in samples rather than seconds and will also return dt in
% samples

% INPUT:
%   basepath            string. path to recording folder {pwd}.
%   dtstr               string of date time if empty will
%                       be extracted from basepath
%   dtFormat            string. format of dtstr. if empty will be extracted
%                       from length of dtstr. best is 'yyMMdd_HHmmss'
%   tstamp              numeric. seconds to add (positive) or remove
%                       (negative) from dtstr record start [seconds]
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
addOptional(p, 'basepath', pwd);
addOptional(p, 'dtstr', []);
addOptional(p, 'dtFormat', []);
addOptional(p, 'tstamp', [], @isnumeric);
addOptional(p, 'fs', 1, @isnumeric);

parse(p, varargin{:})
basepath            = p.Results.basepath;
dtstr               = p.Results.dtstr;
dtFormat            = p.Results.dtFormat;
tstamp              = p.Results.tstamp;
fs                  = p.Results.fs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% conversion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% alternative
% input includes tstamp and a cell of 2 elements, dt. if the first element
% of dt is empty, will be the datetime from basename. if the second element
% is empty or Inf, will be the recording duration. 

% basepath = pwd;
% dtstr = '090000';
% tstamp = [];

% ALT 1: dtstr and tstamp specified
% will add / remove tstamp from dtstr and return output as date time

% ALT 2: dtstr specified but not tstamp 
% find dtstr relative to basename and return in samples / seconds

% if dtstr not specified, will extract it from basename

% if dtstr is Inf, will be the recording duration and tstamp will be taken
% from it

% if fs specified will work in samples rather than seconds 


% if dtformat specified, will try to relate it to either dtstr or basename.
% if not specified will try to extract it from numel of dtstr and basename
% separately (thus both may have different format)

% if dtstr specified as time only (< 7 elements), date will be extracted
% from basename

% find basename date time (assumes mousename_date_time)
[~, basename] = fileparts(basepath);
baseparts = split(basename, '_');
basedate = baseparts{end - 1};
basetime = baseparts{end};
dtstrBase = [basedate '_' basetime];
if numel(dtFormat) == numel(dtstrBase)
    dtBase = datetime(dtstrBase, 'InputFormat', dtFormat);
else
    [dtBase, ~] = guessDateTime(dtstrBase);
end

% find date time of dtstr input
if ~isempty(dtstr)
    if numel(dtstr) < 7
        dtstr = [basedate '_' dtstr];
    end
    if numel(dtFormat) == numel(dtstr)
        dt = datetime(dtstr, 'InputFormat', dtFormat);
    else
        [dt, ~] = guessDateTime(dtstr);
    end
end

% add / remove tstamp from original date time
if tstamp > 0
    dt = dt + seconds(tstamp / fs);
elseif tstamp < 0
    dt = dt - seconds(tstamp / fs);
elseif isempty(tstamp)
    if dt <= dtBase
        dt = dt + hours(24);
    end
    dt = seconds(dt - dtBase) * fs;
end

end

