basepath = 'D:\Data';
basename = 'FsecData.xls';
centerPeak = 10;
win = 2.5;

filename = fullfile(basepath, basename);
raw = readmatrix(filename);
raw(end : -1 : end - 1, :) = [];
ev = raw(:, 2);
fs = 1 / nanmean(diff(ev)); % sampling frequency of elution volume

for i = 3 : size(x, 2)
    
    x = raw(:, i);
    
    % manually find peaks
    % d = diff(x(:, i));
    % idx = find(d(1 : end - 1) >= 0 & d(2 : end) <= 0);
    % idx = idx + 1;
    % plot(x(:, i))
    % hold on
    % scatter(idx, x(idx, i))
    
    % use built-in matlab function
    [p, idx] = findpeaks(x, 'NPeaks', 10, 'MinPeakProminence', 500,...
        'MinPeakDistance', 100);
    [~, idx1] = sort(abs(idx - centerPeak * fs));
    idx = idx(idx1(1 : 2));
        
    % find second peak according to change in slope within a pre-defined
    % window
    % x2 = x(centerPeak * fs : (centerPeak + win) * fs);
    % p = find(ischange(x2, 'linear', 'MaxNumChanges', 5));
%     
%     plot(ev, x)
%     hold on
%     scatter(ev(idx), x(idx))
%     plot(x2)
%     hold on
%     scatter(p, x2(p))
%     
%     [val, posi] = max(x(p, i))
%     
%     plot(diff(d))
    
figure
hold on
plot(x)
yyaxis right
plot(diff(x))
hold on
plot(smooth(diff(diff(x)), 10), '-k')
legend({'x', 'diff', 'diff2'})

end