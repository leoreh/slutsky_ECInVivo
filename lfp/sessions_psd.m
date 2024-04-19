function [bands, psdMat, epochStats] = sessions_psd(varargin)

% arranges the psd in states from different sessions of the same mouse.
% also calculates the power in specific bands. can normalize the psd to
% broadband power and can normalize each session to a baseline session/s.
% plots the results.
%
% INPUT:
%   basepaths       cell of strings describing sessions to analyze. if
%                   empty will use mname
%   mname           string. mouse name to analyze
%   nbins           numeric. divide each recording to nbins {2}
%   idxBsl          numeric. indices to baseline sessions for normalization
%   flgNormTime     logical. normalize each session to baseline (1st
%                   session) {true}. only applied to bands
%   flgNormBand     logical. normalize psd to broadband {true}
%   flgEmg          logical. load [basename].psdEmg.mat {true}
%   flgAnalyze      logical. re-analyze psd {false}
%   flgDb           logical. convert psd to dB {true}
%   saveVar         logical. save ss var {true}
%   saveFig         logical. save figure {true}
%   graphics        logical. plot confusion chart and state separation {true}
%
% OUTPUT:
%   bands           numeric 3d mat of band x state x session
%   powdb           numeric 3d mat of state x freq x session
%
% DEPENDENCIES:
%   calc_psd
%   calc_bands
%
% TO DO LIST:
%
% 20 oct 22 LH      updates
% 31 mar 24             calc_bands and cleanup


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addOptional(p, 'basepaths', {}, @iscell);
addOptional(p, 'mname', '', @ischar);
addOptional(p, 'nbins', 2, @isnumeric);
addOptional(p, 'idxBsl', 1, @isnumeric);
addOptional(p, 'flgNormTime', true, @islogical);
addOptional(p, 'flgNormBand', true, @islogical);
addOptional(p, 'flgEmg', true, @islogical);
addOptional(p, 'flgAnalyze', false, @islogical);
addOptional(p, 'flgDb', false, @islogical);
addOptional(p, 'saveVar', true, @islogical);
addOptional(p, 'saveFig', true, @islogical);
addOptional(p, 'graphics', true, @islogical);

parse(p, varargin{:})
basepaths       = p.Results.basepaths;
mname           = p.Results.mname;
nbins           = p.Results.nbins;
idxBsl          = p.Results.idxBsl;
flgNormTime     = p.Results.flgNormTime;
flgNormBand     = p.Results.flgNormBand;
flgEmg          = p.Results.flgEmg;
flgAnalyze      = p.Results.flgAnalyze;
flgDb           = p.Results.flgDb;
saveVar         = p.Results.saveVar;
saveFig         = p.Results.saveFig;
graphics        = p.Results.graphics;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% vars to load from each session
if flgEmg
    varsFile = ["psdEmg"; "session"];
else
    varsFile = ["psd"; "session"];
end
varsName = ["psd"; "session"];

% load
if isempty(basepaths)
    [v, basepaths] = getSessionVars('mname', mname, 'varsFile', varsFile,...
        'varsName', varsName);
else
    [v, basepaths] = getSessionVars('basepaths', basepaths, 'varsFile', varsFile,...
        'varsName', varsName);
end
nfiles = length(basepaths);

% get mouse name
[~, basename] = fileparts(basepaths{1});
mname = regexp(basename, 'lh\d+', 'match');
mname = mname{1};

% state params
sstates = v(1).psd.info.sstates;
nstates = length(sstates);
snames = v(1).psd.info.snames;
clr = v(1).psd.info.clr;

% psd params
freq = v(1).psd.info.faxis;
spkgrpIdx = v(1).psd.info.spkgrpIdx;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% divide each recordings to nbins and reclaculate vars
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize
psdMat = nan(nfiles * nbins, nstates, length(freq));
timebins = nan(nfiles, nbins, 2);
timepnt = nan(nfiles, 1);               % currently not implemented
cnt = 1;

for ifile = 1 : nfiles
    
    % create timebins
    if isfield(v(ifile).session.general, 'timepnt')
        timepnt(ifile) = v(ifile).session.general.timepnt;
    else
        timepnt(ifile) = nan;
    end
    recLen = floor(v(ifile).session.general.duration);
    timebins(ifile, :, :) = n2chunks('n', recLen, 'nchunks', nbins, 'pnts', []); 
    
    for ibin = 1 : nbins
        for istate = 1 : nstates
            
            % get indices to epochs in timebins
            stateEpochs = v(ifile).psd.epochs.stateEpochs{istate};
            epochIdx = stateEpochs(:, 1) > timebins(ifile, ibin, 1) &...
                stateEpochs(:, 2) < timebins(ifile, ibin, 2);
            
            % recalculate mean psd
            psdMat(cnt, istate, :) =...
                mean(v(ifile).psd.epochs.clean{istate}(epochIdx, :), 1);

            % organize epoch stats
            epLen{istate, cnt} = v(ifile).psd.epochs.epochStats.epLen{istate}(epochIdx);
            nepochs(istate, cnt) = numel(epLen{cnt, istate});
            prctDur(istate, cnt) = sum(epLen{cnt, istate}) * 100 / diff(timebins(ifile, ibin, :));
        end
        cnt = cnt + 1;
    end
end

% convert to dB (chronux output is power not magnitude)
if flgDb
    psdMat = 10 * log10(psdMat);
    ytxt = 'PSD [dB]';
else
    ytxt = 'PSD [mV^2/Hz]';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bands
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% grab band info
[~, info] = calc_bands('psdData', squeeze(psdMat(1, istate, :))',...
    'freq', freq, 'flgNormBand', flgNormBand);
bandNames = info.bandNames;
nbands = length(bandNames);
bandFreqs = info.bandFreqs;

% calc power in specific frequency bands
bands = nan(nfiles * nbins, nstates, nbands);
for istate = 1 : nstates
    [bands(:, istate, :)] = [calc_bands('psdData', squeeze(psdMat(:, istate, :)),...
        'freq', freq, 'flgNormBand', flgNormBand)]';
end

% params if bands are normalized to broadband
if flgNormBand
    sbands = [2 : nbands];
    ytxt_band1 = 'Norm. Power';
else
    sbands = [1 : nbands];
    ytxt_band1 = 'Power [mV^2/Hz]';
end

% normalize bands to baseline
if flgNormTime
    bsl = mean(bands(idxBsl, :, :), 1, 'omitnan');
    bands = (bands ./ bsl) * 100;
    ytxt_band2 = ' (% BSL)';
else
    ytxt_band2 = '';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% organize struct and save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

epochStats.epLen = epLen;
epochStats.nepochs = nepochs;
epochStats.prctDur = prctDur;

if saveVar
    mousepath = fileparts(basepaths{1});
    bandsfile = fullfile(mousepath, [mname, '_psdSessions.mat']);
    save(bandsfile, 'bands')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~graphics
    return
end

% experiment specifics
dlFileIdx = [];         % delete from graphics
idxPer = [4, 11];       % index of perturbation bins
idxPer = [];

% organize colors
clrPsd = zeros(size(psdMat, 1), 3);
if ~isempty(idxPer)
    clrPsd(1 : idxPer(1) - 1, :) = repmat(linspace(0, 0.7, idxPer(1) - 1)', 1, 3);
    clrPsd(idxPer(1) : idxPer(2), :) = autumn(diff(idxPer) + 1);
    clrPsd(idxPer(2) + 1 : end, 3) = linspace(1, 0.7, size(psdMat, 1) - idxPer(2));
else
    clrPsd = turbo(size(psdMat, 1));
end
clrBands = turbo(nbands);

% organize band legend
for iband = 1 : nbands
    lgd{iband} = sprintf('%s [%.1f-%d Hz]',...
        bandNames{iband}, bandFreqs(iband, 1),...
        bandFreqs(iband, 2));
end

% open figure
fh = figure;
setMatlabGraphics(true)
set(fh, 'WindowState','maximized');
tlayout = [5, nstates];
th = tiledlayout(tlayout(1), tlayout(2));
th.TileSpacing = 'tight';
th.Padding = 'none';
title(th, mname, 'interpreter', 'none')

% spectrogram and emg throughout recording
axh = nexttile(th, 1, [1, tlayout(2)]); cla
sessions_catVarTime('dataPreset', 'spec', 'graphics', true,...
    'dataAlt', spkgrpIdx,...
    'basepaths', basepaths, 'xTicksBinsize', 24 / nbins,...
    'markRecTrans', false, 'axh', axh);

% psd per state across time
tilebias = tlayout(2);
for istate = 1 : nstates
    axh = nexttile(th, istate + tilebias, [1, 1]); cla
    hold on
    ydata = squeeze(psdMat(:, istate, :));
    ph = plot(freq, ydata, 'LineWidth', 2);
    set(axh, 'ColorOrder', clrPsd)
    delete(ph(dlFileIdx))
    ph(dlFileIdx) = [];

    set(gca, 'XScale', 'log')
    if ~flgDb
        set(gca, 'YScale', 'log')
    end
    title(snames{istate})
    xlabel('Frequency [Hz]')
    ylabel(ytxt)
    legend(split(num2str(1 : nfiles * nbins)), 'Location', 'Southwest', 'FontSize', 9)
end

% bands across time
tilebias = tlayout(2) * 2;
for istate = 1 : nstates
    axh = nexttile(th, istate + tilebias, [1, 1]); cla
    ph = plot(squeeze(bands(:, istate, sbands)), 'LineWidth', 2);
    set(axh, 'ColorOrder', clrBands(sbands, :))
    hold on
    if flgNormTime
        plot(xlim, [100 100], '--k')
    end
    ylabel([ytxt_band1, ytxt_band2])
end
legend(lgd(sbands))

% epoch stats - epoch length
tilebias = tlayout(2) * 3;
for istate = 1 : nstates
    axh = nexttile(th, istate + tilebias, [1, 1]); cla
    dataMat = cell2nanmat(epLen(:, istate), 2);
    plot_boxMean('dataMat', dataMat, 'clr', clr{istate},...
        'plotType', 'bar', 'axh', axh)
    ylabel('Epoch Length (s)')
end

% epoch stats - percent duration
tilebias = tlayout(2) * 4;
axh = nexttile(th, 1 + tilebias, [1, 1]); cla
ph = plot(prctDur, 'LineWidth', 2);
for istate = 1 : nstates
    ph(istate).Color = clr{istate};
end
ylabel('Duration (%)')
xlim([0.5, nbins * nfiles + 0.5])
xticks([1 : nbins * nfiles])

% epoch stats - percent duration
axh = nexttile(th, 2 + tilebias, [1, 1]); cla
ph = plot(nepochs, 'LineWidth', 2);
for istate = 1 : nstates
    ph(istate).Color = clr{istate};
end
ylabel('No. Epochs')
xlim([0.5, nbins * nfiles + 0.5])
xticks([1 : nbins * nfiles])

% save fig
if saveFig
    figpath = fullfile(mousepath, 'graphics');
    mkdir(figpath)
    figname = fullfile(figpath, [mname, '_psdBands']);
    savefig(fh, figname, 'compact')
end

end

% EOF