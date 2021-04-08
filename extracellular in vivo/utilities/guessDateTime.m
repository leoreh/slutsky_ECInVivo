function [dt, dtFormat] = guessDateTime(dtstr)

% tries to guess the date time format of a string based on its length

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
end
