function [dt, dtFormat] = guessDateTime(dtstr)

% tries to guess the date time format of a string based on its length and
% returns the date time of dtstr. dtstr can also be basename in the format
% "mousename_date_time"

dtstr = char(dtstr);
ln = numel(dtstr);
mflag = 0;
if ln > 13  % assumes basename
    baseparts = split(dtstr, '_');
    basedate = baseparts{end - 1};
    basetime = baseparts{end};
    dtstr = [basedate '_' basetime];
    mflag = 1;
end
ln = numel(dtstr);

if ln == 4
    dtFormat = 'HHmm';
elseif ln == 6
    dtFormat = 'HHmmss';
elseif ln == 11
    dtFormat = 'yyMMdd_HHmm';
elseif ln == 13
    dtFormat = 'yyMMdd_HHmmss';
else
    warning('could not find date time format...\n')
    return
end
dt = datetime(dtstr, 'InputFormat', dtFormat);
if mflag
    dtFormat = ['mname_', dtFormat];
end
end
