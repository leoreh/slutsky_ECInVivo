function [bs, iis, anestates] = aneStates_m(varargin)

% analyze anesthesia states of one mouse. for tetrodes omit the
% basename argument and direct to basepath
%
% INPUT
%   ch          vector. channel/s to analyze.
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
%       getBS
%       getIIS
%       specBand
%       binary2epochs
%
% EXAMPLE
%
% 21 jan 20 LH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addParameter(p, 'ch', 1, @isvector);
addParameter(p, 'basepath', pwd, @isstr);
addParameter(p, 'basename', [], @isstr);
addParameter(p, 'graphics', true, @islogical)
addParameter(p, 'saveVar', true, @islogical);
addParameter(p, 'saveFig', true, @islogical);
addParameter(p, 'forceA', false, @islogical);
addParameter(p, 'forceL', false, @islogical);

parse(p, varargin{:})
ch = p.Results.ch;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lfp data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lfp = getLFP('basepath', basepath, 'ch', ch, 'chavg', {},...
    'fs', 1250, 'interval', [0 inf], 'extension', 'abf', 'pli', true,...
    'savevar', true, 'force', forceL, 'basename', basename);

sig = double(lfp.data(:, ch));
fs = lfp.fs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% burst suppression
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vars = {'std', 'max', 'sum'};
bs = getBS('sig', sig, 'fs', fs, 'basepath', basepath,...
    'graphics', true, 'saveVar', true, 'binsize', 0.5,...
    'clustmet', 'gmm', 'vars', vars, 'basename', basename,...
    'saveFig', false, 'forceA', forceA);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% iis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% thr during suppression
thr(1) = 0;
thr(2) = mean(sig(~bs.binary)) + 15 * std(sig(~bs.binary));

marg = 0.1;
iis = getIIS('sig', sig, 'fs', fs, 'basepath', basepath,...
    'graphics', true, 'saveVar', saveVar, 'binsize', 300,...
    'marg', marg, 'basename', basename, 'thr', thr,...
    'saveFig', false, 'forceA', forceA);
wvstamps = -marg * 1000 : 1 / fs * 1000 : marg * 1000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% burst suppression after IIS filteration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield(iis, 'filtered') && size(iis.wv, 1) > 50
    bs = getBS('sig', iis.filtered, 'fs', fs, 'basepath', basepath,...
        'graphics', false, 'saveVar', true, 'binsize', 1,...
        'clustmet', 'gmm', 'vars', vars, 'basename', basename,...
        'saveFig', false, 'forceA', forceA);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spectrum
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calc dband
[dband, tband] = specBand('sig', sig, 'graphics', false, 'band', [1 4]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IIS in anesthesia states
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
anestates = [];
% find epochs
% limbs = [0 0.5];
% limd = [0.5 1];
% vec = bs.bsr' > limbs(1) & bs.bsr' < limbs(2) &...
%     dband' > limd(1) & dband' < limd(2);
% 
% % binsize of bsr and delta power is 10 s. minimum duration is 1 min [6
% % bins] and minimum time between epochs of deep anesthesia is 1 minutes
% % [6 bins]
% epochs = binary2epochs('vec', vec, 'minDur', 6, 'interDur', 6);
% 
% % epoch bins to sample idx
% epochs = bs.cents(epochs);
% 
% % find IIS within anesthesia epochs
% idx = [];
% for j = 1 : size(epochs, 1)
%     idx = [idx; find(iis.peakPos > epochs(j, 1) &...
%         iis.peakPos < epochs(j, 2))];
% end
% wv = iis.wv(idx, :);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if graphics
    
    f = figure;
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    suptitle(basename)
    
    % raw
    sb1 = subplot(3, 2, 1 : 2);
    plot((1 : length(sig)) / fs / 60, sig, 'k')
    axis tight
    ylim([-1 2])
    hold on
    plot(xlim, [iis.thr(2) iis.thr(2)], '--r')
    ylabel('Voltage [mV]')
    yyaxis right
    p1 = plot(iis.cents / fs / 60, iis.rate, 'r', 'LineWidth', 3);
    if ~isempty(p1)
        p1.Color(4) = 0.3;
    end
    ylim([0 1])
    ylabel('Rate [IIS / min]')
    legend({'Raw', 'IIS thr', 'IIS rate'}, 'Location', 'northwest')
    axis tight
    set(gca, 'TickLength', [0 0])
    box off
    
    % bsr and delta
    sb2 = subplot(3, 2, 3 : 4);
    plot(bs.cents / fs / 60, bs.bsr, 'k', 'LineWidth', 2)
    hold on
    plot(tband / 60, dband, 'b', 'LineWidth', 2)
    legend({'BSR', 'Delta'})
    Y = ylim;
    if exist('epochs', 'var')
        fill([epochs fliplr(epochs)]' / fs / 60, [Y(1) Y(1) Y(2) Y(2)],...
            'b', 'FaceAlpha', 0.2,  'EdgeAlpha', 0, 'HandleVisibility', 'off');
    end
    xlabel('Time [m]')
    ylabel('[a.u.]')
    axis tight
    set(gca, 'TickLength', [0 0])
    box off
    title('Surgical Anesthesia')
    
    % bs
    subplot(3, 2, 5);
    midsig = round(length(sig) / 2 / fs / 60);
    ylimsig = sig(midsig * fs * 60 : (midsig + 2) * fs * 60);
    plot((1 : length(sig)) / fs / 60, sig)
    hold on
    axis tight
    ylim([min(ylimsig) max(ylimsig)])
    Y = ylim;
    fill([bs.stamps fliplr(bs.stamps)] / fs / 60, [Y(1) Y(1) Y(2) Y(2)],...
        'k', 'FaceAlpha', 0.4,  'EdgeAlpha', 0);
    xlim([midsig midsig + 2])
    ylabel('Voltage [mV]')
    set(gca, 'TickLength', [0 0])
    box off
    title('BS')
    
    % iis waveforms
    if ~isempty(iis.wv)
        subplot(3, 2, 6)
        plot(wvstamps, iis.wv)
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

