function [iis] = getIIS(varargin)

% detects inter-ictal spikes from LFP.
% 
% INPUT
%   sig         signal for detection
%   fs          sampling frequency
%   binsize     scalar {60} in [s]. for rate calculation
%   marg        scalar {0.1} in [s]. time margin for clipping spikes 
%   thr         scalar {15} in [z scores]. used for detection 
%   basepath    recording session path {pwd}
%   basename    string. if empty extracted from basepath
%   graphics    logical {true}. plot figure
%   saveVar     logical {true}. save variable
%   saveFig     logical {true}. save figure
%   forceA      logical {false}. force analysis even if .mat exists
%
% OUTPUT
%   iis         struct with fields:
%       stamps      n x 2 mat where n is the number of events and the
%                   columns represent the start and end of each event [samps]
%       peaks       the timestamp for each peak [samps]
%       peakPower   the voltage at each peak [uV]
%       <params>    as in input
%
% TO DO LIST
%       # investigate thr
%       # filter IIS by removing average waveform
%       # find width via wavelet and after hyperpolirazation
%
% CALLS
%       calcFR     
%
% 02 jan 20 LH.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addParameter(p, 'sig', [], @isnumeric)
addParameter(p, 'fs', 1250, @isnumeric)
addParameter(p, 'binsize', 300, @isnumeric)
addParameter(p, 'marg', [], @isnumeric)
addParameter(p, 'thr', 15, @isnumeric)
addParameter(p, 'basepath', pwd, @isstr);
addParameter(p, 'basename', [], @isstr);
addParameter(p, 'graphics', true, @islogical)
addParameter(p, 'saveVar', true, @islogical);
addParameter(p, 'saveFig', true, @islogical);
addParameter(p, 'forceA', false, @islogical);

parse(p, varargin{:})
sig = p.Results.sig;
fs = p.Results.fs;
marg = p.Results.marg;
thr = p.Results.thr;
binsize = p.Results.binsize;
basepath = p.Results.basepath;
basename = p.Results.basename;
graphics = p.Results.graphics;
saveVar = p.Results.saveVar;
saveFig = p.Results.saveFig;
forceA = p.Results.forceA;

if isempty(marg)
    marg = round(0.1 * fs);
end

wvstamps = -marg / fs * 1000 : 1 / fs * 1000 : marg / fs * 1000;
tstamps = [1 : length(sig)] / fs;

binsize = binsize * fs;     % because sig to calcFR is in samples

% initialize output
iis = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check if mat already exists
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(basename)
    [~, basename] = fileparts(basepath);
end
filename = [basepath, '\', basename, '.iis.mat'];
if ~forceA
    if exist(filename)
        load(filename)
        fprintf('\n loading %s \n\n', filename)
        return
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prep signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = sig * 1000;         % so most values are above 1
x = zscore((x / prctile(x, 99)) .^ 2);

% investigate threshold
% j = 1;
% idx = 5 : 0.5 : floor(max(x));
% nevents = zeros(1, length(idx));
% for i = 1 : length(idx)
%     nevents(i) = sum(x > idx(i));
% end
% figure
% plot(idx, cumsum(nevents))
% yyaxis right
% plot(idx, cumsum(log10(nevents)))
% axis tight
% 
% histogram((nevents), 50, 'Normalization', 'cdf')

        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% detect spikes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iie = find([0 diff(x > thr(1)) > 0] & abs(sig) > 1);

% return if no IIS are detected
nspks = length(iie);
if isempty(iie)
    disp('\nno inter-ictal spikes found\n');
    return
else
    disp(sprintf('\n%d inter-ictal spikes \n', nspks));
end

% thr in mV
[~, idx] = min(abs(x - thr(1)));
thr(2) = sig(idx);

% select local maximum and clip spikes
peak = zeros(length(iie), 1);
pos = zeros(length(iie), 1);
seg = zeros(length(iie), marg * 2 + 1);
for i = 1 : length(iie)
    seg(i, :) = sig(iie(i) - marg : iie(i) + marg);
    [peak(i), pos(i)] = max(abs(seg(i, :)));
    pos(i) = iie(i) - marg + pos(i);
    seg(i, :) = sig(pos(i) - marg : pos(i) + marg);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% power spectrum via wavelet
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% filter bank
fb = cwtfilterbank('SignalLength', size(seg, 2), 'VoicesPerOctave', 32,...
    'SamplingFrequency', fs, 'FrequencyLimits', [1 625]);

% average wavelet coefficients for each spike. note this produces very
% similar results to the coefficients obtained from the average spike
% waveform but is still carried out since the wavelet of each spike is
% already calculated for measureing width
for i = 1 : size(seg, 1)
    [cfs(i, :, :), f, coi] = cwt(seg(i, :), 'FilterBank', fb);

    % spike width by inverse of max wavelet
    [mp, ifreq] = max(abs(squeeze(cfs(i, :, :))), [], 1);
    [~, ip] = max(mp);
    ifreq = ifreq(ip);
    iis.maxf(i) = f(ifreq);
    iis.spkw(i) = 1000 / iis.maxf(i);
    
    % spike width at half maximum
    % [~, ~, w(i)] = findpeaks(-seg(i, :), fs, 'NPeaks', 1,...
    %    'WidthReference', 'halfheight', 'MinPeakHeight', abs(thr(2)));    
end
cfs = squeeze(mean(abs(cfs), 1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rate 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[iis.rate, iis.edges, iis.cents] = calcFR(pos, 'winCalc', [1, length(sig)],...
    'binsize', binsize, 'smet', 'MA', 'c2r', false);
iis.rate = iis.rate / fs * 60;      % convert to minutes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ACG 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ccg, tccg] = CCG({pos / fs}, [], 'duration', 10,...
    'binSize', 0.1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arrange struct output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iis.peakPos = pos;
iis.peakPower = peak;
iis.fs = fs;
iis.binsize = binsize;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if saveVar
    save(filename, 'iis')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if graphics
    fh = figure;
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    suptitle(basename)
    
    % raw and iis
    subplot(3, 4, 1 : 2)
    plot([pos pos] / fs / 60, [min(sig), max(sig)], '--g', 'LineWidth', 0.1)
    hold on
    plot(tstamps / 60, sig)
    hold on
    axis tight
    ylabel('Voltage [mV]')
    set(gca, 'TickLength', [0 0])
    box off
    title('Raw signal')
    
    % detection
    subplot(3, 4, 5 : 6)
    plot(tstamps / 60, x, 'k')
    hold on
    axis tight
    plot(xlim, [thr(1) thr(1)], '--r')
    plot([pos pos] / fs / 60, [-10 -1], '--g', 'LineWidth', 2)
    ylabel('Z-score')
    set(gca, 'TickLength', [0 0])
    box off
    title('IIS detection')
    
    % iis rate
    subplot(3, 4, 9 : 10)
    plot(iis.cents / fs / 60, iis.rate, 'k', 'LineWidth', 3)
    ylabel('Rate [spikes / min]')
    xlabel('Time [m]')
    axis tight
    set(gca, 'TickLength', [0 0])
    box off
    title('IIS rate')
    
    % iis waveforms
    subplot(3, 4, 3)
    plot(wvstamps, seg')
    ylabel('Voltage [mV]')
    xlabel('Time [ms]')
    axis tight
    xticks([-marg, 0, marg] / fs * 1000);
    set(gca, 'TickLength', [0 0])
    box off
    title('Spike waveform')
    
    % mean + std waveform
    axes('Position',[.542 .71 .09 .07])
    box on
    stdshade(seg, 0.5, 'k', wvstamps)
    axis tight
    set(gca, 'TickLength', [0 0], 'YTickLabel', [], 'XTickLabel', [],...
        'XColor', 'none', 'YColor', 'none', 'Color', 'none')
    box off
    
    % iis cwt
    subplot(3, 4, 4)
    imagesc(wvstamps, f, abs(cfs))
    axis xy
    set(gca, 'YScale', 'log')
    % caxis([0.0001 1])
    origSize = get(gca, 'Position');
    colorbar
    set(gca, 'Position', origSize);
    hold on
    plot(wvstamps, coi, 'w', 'LineWidth', 2)
    ylim([min(f) max(f)]);
    xlim([wvstamps(1) wvstamps(end)]);
    yticks([ceil(min(f)), 50, 100, max(f)])
    xticks([min(wvstamps), 0, max(wvstamps)])
    ylabel('Frequency [Hz]')
    xlabel('Time [ms]')
    set(gca, 'TickLength', [0 0])
    box off
    title('Average Scalogram')
      
    % amplitude histogram
    subplot(3, 4, 7)
    h = histogram(log10(iis.spkw), 20, 'Normalization', 'Probability');
    h.EdgeColor = 'none';
    h.FaceColor = 'k';
    h.FaceAlpha = 1;
    xlabel('Spike Width [ms]')
    ylabel('Probability [%]')
    set(gca, 'TickLength', [0 0])
    box off
    title('IIS width')
    
    % amplitude histogram
    subplot(3, 4, 8)
    h = histogram(log10(iis.peakPower), 20, 'Normalization', 'Probability');
    h.EdgeColor = 'none';
    h.FaceColor = 'k';
    h.FaceAlpha = 1;
    xlabel('Peak voltage [log(uV)]')
    ylabel('Probability [%]')
    set(gca, 'TickLength', [0 0])
    box off
    title('IIS amplitude')
        
    % max frequency and amplitude vs. time
    subplot(3, 4, 11)
    s1 = scatter(iis.peakPos / fs / 60, iis.maxf, 2, 'k', 'filled');
    l1 = lsline;
    ylabel('Frequency [Hz]')
    set(l1, 'color', 'k')
    l1.LineWidth = 3;
    yyaxis right
    s2 = scatter(iis.peakPos / fs / 60, iis.peakPower, 2, 'b', 'filled');
    l2 = lsline;
    ylabel('Amplitude [mV]')
    set(l2, 'color', 'b')
    l2.LineWidth = 3;
    yyaxis right
    axis tight
    xlabel('Time [m]')
    set(gca, 'TickLength', [0 0])
    box off
    title('Frequency. vs. Time')
    
    % amplitude vs. time
    subplot(3, 4, 12)
    plotCCG('ccg', ccg, 't', tccg / 1000, 'basepath', basepath,...
        'saveFig', false, 'c', {'k'});
    xlabel('Time [s]')
    title('Autocorrelogram')
    
    if saveFig
        figname = [basename '_IIS'];
        export_fig(figname, '-tif', '-transparent')
        % savePdf(figname, basepath, ff)
    end
end
end

% EOF

