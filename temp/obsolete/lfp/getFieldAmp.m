function [amp, rm] = getFieldAmp(varargin)

% loads data specified in filename. Each file should contain traces with
% the same stimulus intensity. allows user to remove unwanted traces.
% plots amplitude as a function of stimulus intensity and superimposed
% traces.
%
% INPUT
%   sig         signal mat. voltage (rows) x traces (columns)
%   fs          sampling frequency
%   nstim       number of stimulations in one trace.
%   start       start times of stimulations [s]. if empty extracted from
%               artifacts. must be equal in length to nstim.
%   stop        times of maximum response [s]. if empty extracted from min
%               within window. must be equal in length to nstim
%   inspect     logical. inspect traces {1} or not (0).
%   basepath    recording session path {pwd} to save figure and variables
%   filename    filename to save fig and var. if empty that datetime
%   graphics    logical. plot graphics {1} or not.
%   saveFig     logical. saveFig to current path {1} or not (0).
%   saveVar     logical. save output to current path {1} or not (0).
%
% OUTPUT
%   amp         mat of amp responses for each stim (rows) x trace (column)
%   rm          index of traces that were removed after inspection
%
% TO DO LIST
%   # find start times from max(diff) or corr
%
% 12 feb 20 LH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addOptional(p, 'sig', []);
addOptional(p, 'fs', 100300.902708124, @isnumeric);
addOptional(p, 'start', [], @isnumeric);
addOptional(p, 'stop', [], @isnumeric);
addOptional(p, 'inspect', true, @islogical);
addOptional(p, 'nstim', 5);
addOptional(p, 'basepath', pwd, @isstr);
addOptional(p, 'filename', datestr(datetime), @isstr);
addOptional(p, 'graphics', true, @islogical);
addOptional(p, 'saveFig', false, @islogical);
addOptional(p, 'saveVar', false, @islogical);

parse(p,varargin{:})
sig = p.Results.sig;
fs = p.Results.fs;
start = p.Results.start;
stop = p.Results.stop;
inspect = p.Results.inspect;
nstim = p.Results.nstim;
basepath = p.Results.basepath;
filename = p.Results.filename;
graphics = p.Results.graphics;
saveFig = p.Results.saveFig;
saveVar = p.Results.saveVar;

% constants
marg = round(0.002 * fs);
t2res = round(0.015 * fs);
t2max = round(0.005 * fs);
tstamps = [1 : size(sig, 1)] / fs;

% output
rm = 0;
amp = nan(nstim, size(sig, 2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prepare data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

avgsig = mean(sig, 2);

% remove DC calculated on mean of all traces concatenated
sz = size(sig);
sig = rmDC(sig(:));
sig = reshape(sig, sz);

% manually inspect and remove unwanted traces
if inspect
    [sig, rm] = rmTraces(sig, tstamps);
end

% find stim start times [in samples]
if isempty(start)
    % ALT 1 - use pairwise correlations between traces
    %     artcor = getEMGfromLFP(sig, 'emgFs', fs / 100, 'fs', fs,...
    %         'saveVar', false, 'graphics', false);
    %     [~, idx] = sort(artcor.data, 'descend');
    %     x = idx(1 : 5);
    %     x = x(:) - x(:)';
    %     diag(x)
    % ALT 2 - use max(diff).
else
    % find start from max value in win after artifact
    start = round(start * fs + marg);
    for i = 1 : nstim
        win = start(i) : start(i) + t2max;
        [~, idx(i)] = max(avgsig(win));
        start(i) = start(i) + idx(i);
    end
end

% find stop times [in samples] according to min value in win after artifact
if isempty(stop)
    for i = 1 : nstim
        win = start(i) + marg : round(start(i) + t2res);
%         amp(i, :) = abs(min(sig(win, :))) - abs(sig(start(i), :));
        [~, idx(i)] = min(avgsig(win, :));
        stop(i) = win(1) + idx(i);
    end
else
    stop = stop * fs;
end

% find amp response according to diff between start and stop
for i = 1 : nstim
    amp(i, :) = sig(start(i), :) - (sig(stop(i), :));
end


if saveVar
    save(filename, 'amp')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if graphics
   
    fh = figure;
    subplot(2, 2, 1)
    plot(tstamps, sig)
    ylabel('Voltage [mV]')
    set(gca, 'TickLength', [0 0])
    box off
    axis tight
    title('Traces')
   
    subplot(2, 2, 3)
    stdshade(sig', 0.5, 'k', tstamps)
    ylabel('Voltage [mV]')
    xlabel('Time [s]')
    set(gca, 'TickLength', [0 0])
    box off
    axis tight
   
    subplot(2, 2, 2)  
    boxplot(amp', 'BoxStyle', 'outline', 'Color', 'k', 'notch', 'off')
    hold on
    for i = 1 : nstim
        scatter(ones(1, size(sig, 2)) * i, amp(i, :), 20, 'b', 'filled')
    end
    ylabel('Amplitude [mV]')
    xticks([1 : nstim])
    set(gca, 'TickLength', [0 0])
    box off
    title('Amplitude')
   
    subplot(2, 2, 4)
    plot(mean(amp, 2) / mean(amp(1, :)))
    ylabel('Norm. Amplitude')
    xlabel('Stim #')
    xticks([1 : nstim])
    set(gca, 'TickLength', [0 0])
    box off
   
    if saveFig
        savePdf(filename, basepath, fh)
    end
end

end