function [bs, iis, ep] = aneStates_m(varargin)

% analyze anesthesia states of one mouse. for tetrodes omit the
% basename argument and direct to basepath
%
% INPUT
%   ch          vector. channel/s to analyze.
%   binsize     scalar [s]. applied to bsr, iis rate, and spec band
%   smf         smooth factor [bins]. applied to bsr, iis rate, and spec band
%   basepath    string. recording session path {pwd}
%   basename    string. mouse basename. if empty extracted from basepath
%   graphics    logical. plot figure {1}
%   saveVar     logical. save variable {1}
%   saveFig     logical. save figure {1}
%   forceA      logical {0}. force analysis even if .mat exists
%   forceL      logical {0}. force load lfp even if .mat exists
%
% OUTPUT
%   bs          struct (see getBS.m)
%   iis         struct (see getIIS.m)
%   spec        struct with fields:
%       dband   normalized power in the delta band (1-4 Hz)
%       sband   normalized power in the sigma band (12-20 Hz)
%       tband   timestamps for normalized power
%
% TO DO LIST
%       #
%
% CALLS
%       fetLFP
%       getBS
%       getIIS
%       specBand
%       binary2epochs
%
% EXAMPLE
%       [bs, iis, ep] = aneStates_m('ch', 1, 'basepath', basepath,...
%        'basename', basename, 'graphics', graphics, 'saveVar', saveVar,...
%        'saveFig', saveFig, 'forceA', forceA, 'binsize', 30, 'smf', 6);
% 
% 21 jan 20 LH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addParameter(p, 'ch', 1, @isvector);
addParameter(p, 'binsize', 60, @isnumeric);
addParameter(p, 'smf', 6, @isnumeric);
addParameter(p, 'basepath', pwd, @isstr);
addParameter(p, 'basename', [], @isstr);
addParameter(p, 'graphics', true, @islogical)
addParameter(p, 'saveVar', true, @islogical);
addParameter(p, 'saveFig', true, @islogical);
addParameter(p, 'forceA', false, @islogical);
addParameter(p, 'forceL', false, @islogical);

parse(p, varargin{:})
ch = p.Results.ch;
binsize = p.Results.binsize;
smf = p.Results.smf;
basepath = p.Results.basepath;
basename = p.Results.basename;
graphics = p.Results.graphics;
saveVar = p.Results.saveVar;
saveFig = p.Results.saveFig;
forceA = p.Results.forceA;
forceL = p.Results.forceL;

if isempty(basename)
    [~, basename] = fileparts(basepath);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lfp data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lfp = getLFP('basepath', basepath, 'ch', ch, 'chavg', {},...
    'fs', 1250, 'interval', [0 inf], 'extension', 'abf', 'pli', true,...
    'savevar', true, 'force', forceL, 'basename', basename);
sig = double(lfp.data(:, ch));
fs = lfp.fs;
binsize = (2 ^ nextpow2(binsize * fs));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% burst suppression
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vars = {'std', 'max', 'sum'};
bs = getBS('sig', sig, 'fs', fs, 'basepath', basepath, 'graphics', graphics,...
    'saveVar', saveVar, 'binsize', 0.5, 'BSRbinsize', binsize, 'smf', smf,...
    'clustmet', 'gmm', 'vars', vars, 'basename', basename,...
    'saveFig', saveFig, 'forceA', forceA, 'vis', false);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% iis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% thr during suppression
thr = [0 mean(sig(~bs.binary(1 : 15 * fs * 60)))...
    + 15 * std(sig(~bs.binary(1 : 15 * fs * 60)))];
marg = 0.05;
iis = getIIS('sig', sig, 'fs', fs, 'basepath', basepath,...
    'graphics', graphics, 'saveVar', saveVar, 'binsize', binsize,...
    'marg', marg, 'basename', basename, 'thr', thr, 'smf', 7,...
    'saveFig', saveFig, 'forceA', forceA, 'spkw', false, 'vis', false);
wvstamps = linspace(-marg, marg, floor(marg * fs) * 2 + 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% burst suppression after IIS filteration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield(iis, 'filtered') && size(iis.wv, 1) > 50
    bs = getBS('sig', iis.filtered, 'fs', fs, 'basepath', basepath, 'graphics', graphics,...
        'saveVar', saveVar, 'binsize', 0.5, 'BSRbinsize', binsize, 'smf', smf,...
        'clustmet', 'gmm', 'vars', vars, 'basename', basename,...
        'saveFig', saveFig, 'forceA', forceA, 'vis', false);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spectrum power bands
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[ep.dband, tband] = specBand('sig', sig, 'graphics', false,...
    'band', [1 4], 'binsize', binsize, 'smf', smf);
[ep.sband, ~] = specBand('sig', sig, 'graphics', false,...
    'band', [12 20], 'binsize', binsize, 'smf', smf);
[ep.gband, ~] = specBand('sig', sig, 'graphics', false,...
    'band', [30 80], 'binsize', binsize, 'smf', smf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IIS in anesthesia states
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% episodes of BS
vec = [bs.bsr > 0.3 & bs.bsr < 0.8 & ep.dband > 0.5];
ep.bs_stamps = binary2epochs('vec', vec, 'minDur', 1, 'interDur', 1);

% bins to samples
ep.bs_stamps = bs.cents(ep.bs_stamps);
if ~isempty(ep.bs_stamps)
    if ep.bs_stamps(1) == bs.cents(1)
        ep.bs_stamps(1) = 1;
    end
    if ep.bs_stamps(end) == bs.cents(end)
        ep.bs_stamps(end) = length(sig);
    end
end

idx = [];
idx2 = [];
for j = 1 : size(ep.bs_stamps, 1)
    % iis within epochs
    idx = [idx; find(iis.peakPos > ep.bs_stamps(j, 1) &...
        iis.peakPos < ep.bs_stamps(j, 2))];
    % mean delta within epochs
    idx2 = [idx2, find(iis.cents >= ep.bs_stamps(j, 1) &...
        iis.cents <= ep.bs_stamps(j, 2))];
end
ep.bs_nspks = length(idx);
ep.bs_delta = mean(ep.dband(idx2));
wv = iis.wv(idx, :);

% episodes of B
ep.b_stamps = binary2epochs('vec', bs.bsr < 0.3,...
    'minDur', 1, 'interDur', 1);
ep.b_stamps = bs.cents(ep.b_stamps);
if ~isempty(ep.b_stamps)
    if ep.b_stamps(1) == bs.cents(1)
        ep.b_stamps(1) = 1;
    end
    if ep.b_stamps(end) == bs.cents(end)
        ep.b_stamps(end) = length(sig);
    end
end
idx = [];
idx2 = [];
for j = 1 : size(ep.b_stamps, 1)
    idx = [idx; find(iis.peakPos > ep.b_stamps(j, 1) &...
        iis.peakPos < ep.b_stamps(j, 2))];   
    idx2 = [idx2, find(iis.cents >= ep.b_stamps(j, 1) &...
        iis.cents <= ep.b_stamps(j, 2))];
end
ep.b_delta = mean(ep.dband(idx2));
ep.b_nspks = length(idx);

% episode and recording duration
ep.bDur = sum(diff(ep.b_stamps, [], 2));
ep.bsDur = sum(diff(ep.bs_stamps, [], 2));
ep.recDur = length(sig);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if graphics    
    fh = figure('Visible', 'on');
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    suptitle(basename)
    tstamps = (1 : length(sig)) / fs / 60;
    
    % raw
    sb1 = subplot(3, 2, 1 : 2);
    plot(tstamps, sig, 'k')
    axis tight
%     ylim([-1 2])
    hold on
    plot(xlim, [iis.thr(2) iis.thr(2)], '--r')
    ylabel('Voltage [mV]')
    yyaxis right
    p1 = plot(iis.cents / fs / 60, iis.rate, 'r', 'LineWidth', 3);
    if ~isempty(p1)
        p1.Color(4) = 0.3;
    end
    ylim([0 1])
    ylabel('Rate [spikes / bin]')
    legend({'Raw', 'IIS thr', 'IIS rate'}, 'Location', 'northwest')
    axis tight
    set(gca, 'TickLength', [0 0])
    box off
    
    % bsr and delta
    sb2 = subplot(3, 2, 3 : 4);
    plot(bs.cents / fs / 60, bs.bsr, 'k', 'LineWidth', 2)
    hold on
    plot(tband / 60, ep.dband, 'b', 'LineWidth', 2)
    plot(tband / 60, ep.sband, 'r', 'LineWidth', 2)
    plot(tband / 60, ep.gband, 'g', 'LineWidth', 2)
    legend({'BSR', 'Delta', 'Sigma', 'Gamma'})
    Y = ylim;
    if ~isempty(ep.bs_stamps)
        fill([ep.bs_stamps fliplr(ep.bs_stamps)]' / fs / 60, [Y(1) Y(1) Y(2) Y(2)],...
            'b', 'FaceAlpha', 0.2,  'EdgeAlpha', 0, 'HandleVisibility', 'off');
    end
    if ~isempty(ep.b_stamps)
        fill([ep.b_stamps fliplr(ep.b_stamps)]' / fs / 60, [Y(1) Y(1) Y(2) Y(2)],...
            'g', 'FaceAlpha', 0.2,  'EdgeAlpha', 0, 'HandleVisibility', 'off');
    end
    xlabel('Time [m]')
    ylabel('[a.u.]')
    axis tight
    set(gca, 'TickLength', [0 0])
    box off
    title('Anesthesia States')
    
    % zoom in raw, bs, and iis
    subplot(3, 2, 5);
    minmarg = 1.5;
    midsig = round(length(sig) / 2);
    idx = round(midsig - minmarg * fs * 60 : midsig + minmarg * fs * 60);
    idx2 = iis.peakPos > idx(1) & iis.peakPos < idx(end);
    plot(tstamps(idx), sig(idx), 'k')
    axis tight
    hold on
    x = xlim;
    plot(x, [iis.thr(2) iis.thr(2)], '--r')
    scatter(iis.peakPos(idx2) / fs / 60,...
        iis.peakPower(idx2), '*');
    bsstamps = RestrictInts(bs.stamps, [idx(1) idx(end)]);
    Y = ylim;
    if ~isempty(bsstamps)
        fill([bsstamps fliplr(bsstamps)] / fs / 60, [Y(1) Y(1) Y(2) Y(2)],...
            'k', 'FaceAlpha', 0.25,  'EdgeAlpha', 0);
    end
    ylabel('Voltage [mV]')
    xlabel('Time [m]')
    xticks([ceil(midsig / fs / 60 - minmarg), floor(midsig / fs / 60 + minmarg)])
    set(gca, 'TickLength', [0 0])
    box off
    title('IIS')    
    
    % iis waveforms
    if ~isempty(iis.wv)
        subplot(3, 2, 6)
        plot(wvstamps * 1000, iis.wv)
        ylabel('Voltage [mV]')
        xlabel('Time [ms]')
        axis tight
        xticks([-marg, 0, marg] * 1000);
        set(gca, 'TickLength', [0 0])
        box off
        title('IIS waveform')       
        
        % mean + std waveform
        axes('Position',[.571 .11 .15 .1])
        box on
        stdshade(iis.wv, 0.5, 'k', wvstamps)
        axis tight
        set(gca, 'TickLength', [0 0], 'YTickLabel', [], 'XTickLabel', [],...
            'XColor', 'none', 'YColor', 'none', 'Color', 'none')
        title(sprintf('n = %d', size(iis.wv, 1)));
        box off
    end
    
    linkaxes([sb1, sb2], 'x');
    
    if saveFig
        figname = [basename];
        export_fig(figname, '-tif', '-transparent')
        % savePdf(figname, basepath, ff)
    end  
end
end

