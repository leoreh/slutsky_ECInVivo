function psd = psd_states(varargin)

% wrapper for calc_psd dedicated to sleep states. uses the eeg signal from
% sSig by default, but can load any specified ch from a binary file.
% if sleep_states.mat is not found, will separate the recording to "AW" and
% "NREM" according to high- and low-emg activity.
%
% INPUT:
%   basepath        string. path to recording folder {pwd}
%   sig             vector of lfp signal. if empty will grab from sSig or
%                   sigfile (see below)
%   ch              numeric. channels to load from the sigfile. if a vector
%                   is specified, the signal will be averaged across
%                   channels
%   nchans          numeric. no. channels in sigfile
%   fs              numeric. sampling frequency of lfp data. if empty will
%                   be extracted from session struct. also determines the
%                   new sampling frequency of the eeg signal.
%   sigfile         char. name of file to load signal from. if empty but
%                   channel is specified, will load from [basename.lfp]
%   btimes          cell of n x 2 mats. if empty will load from sleep states
%   sstates         numeric. index of selected states to calculate psd
%   ftarget         numeric. requested frequencies for calculating the psd
%   flgEmg          logical. calc psd in high- and low-emg even if states
%                   file exists
%   flgPli          logical. interpret freqs of power line interference (true)
%   saveVar         logical. save ss var {true}
%   forceA          logical. reanalyze recordings even if ss struct
%                   exists (false)
%   graphics        logical. plot confusion chart and state separation {true}
%
% DEPENDENCIES:
%   calc_psd
%
% TO DO LIST:
%
% 07 sep 22 LH  updates:
% 20 mar 24 LH      cleanup, emg params, band analysis
% 26 mar 24 LH      removed win functionality since psd is now calculated
%                   per bout. much simpler this way. bkup exists in temp
%                   folder
% 01 feb 25 LH      added pli interpolation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addOptional(p, 'basepath', pwd);
addOptional(p, 'sig', [], @isnumeric);
addOptional(p, 'ch', [], @isnumeric);
addOptional(p, 'nchans', [], @isnumeric);
addOptional(p, 'fs', 1, @isnumeric);
addOptional(p, 'sigfile', [], @ischar);
addOptional(p, 'sstates', [1, 4, 5], @isnumeric);
addOptional(p, 'btimes', []);
addOptional(p, 'ftarget', 0.5 : 0.5 : 100, @isnumeric);
addOptional(p, 'flgEmg', false, @islogical);
addOptional(p, 'flgPli', true, @islogical);
addOptional(p, 'saveVar', true, @islogical);
addOptional(p, 'forceA', false, @islogical);
addOptional(p, 'graphics', true, @islogical);

parse(p, varargin{:})
basepath        = p.Results.basepath;
sig             = p.Results.sig;
ch              = p.Results.ch;
nchans          = p.Results.nchans;
fs              = p.Results.fs;
sigfile         = p.Results.sigfile;
sstates         = p.Results.sstates;
btimes          = p.Results.btimes;
ftarget         = p.Results.ftarget;
flgEmg          = p.Results.flgEmg;
flgPli          = p.Results.flgPli;
saveVar         = p.Results.saveVar;
forceA          = p.Results.forceA;
graphics        = p.Results.graphics;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% graphics flags
flgSaveFig = true;

% file
cd(basepath)
[~, basename] = fileparts(basepath);
sleepfile = fullfile(basepath, [basename, '.sleep_sig.mat']);

% load session vars
varsFile = ["sleep_states"; "datInfo"; "session"];
varsName = ["ss"; "datInfo"; "session"];
v = getSessionVars('basepaths', {basepath}, 'varsFile', varsFile,...
    'varsName', varsName);

% state params
if flgEmg
    namePrefix = 'psdEmg';
    sstates = [1 : 2];
    cfg = as_loadConfig('flgEmg', true);
    statesfile = fullfile(basepath, [basename, '.sleep_statesEmg.mat']);
else
    namePrefix = 'psd';
    cfg = as_loadConfig();
    if isempty(sstates)
        sstates = 1 : cfg.nstates;
    end
    statesfile = fullfile(basepath, [basename, '.sleep_states.mat']);
end

clr = cfg.colors(sstates);
snames = cfg.names(sstates);
nstates = length(sstates);

% check if already analyzed
psdfile = fullfile(basepath, [basename, '.', namePrefix, '.mat']);
if exist(psdfile, 'file') && ~forceA
    load(psdfile)
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load signal from sSig or from binary if ch specified
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(sig)
    if ~isempty(ch)         % from binary

        % if no file specified, load from .lfp binary
        if isempty(sigfile)
            sigfile = fullfile(basepath, [basename, '.lfp']);
            fs = v.session.extracellular.srLfp;
            nchans = v.session.extracellular.nChannels;
        end

        sig = double(bz_LoadBinary(sigfile,...
            'duration', Inf,...
            'frequency', fs, 'nchannels', nchans, 'start', 0,...
            'channels', ch, 'downsample', 1));

        % average tetrode
        if length(ch) > 1
            sig = mean(sig, 2);
        end

    else                    % from sSig

        sig = load(sleepfile, 'eeg');
        sig = sig.eeg;
        load(sleepfile, 'fs');
        load(sleepfile, 'info');
        ch = info.eegCh;
    end
else                        % assumes input signal is from sSig
    load(sleepfile, 'info');
    ch = info.eegCh;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get state bouts if not provided
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(btimes)
    ss = load(statesfile);
    varName = fieldnames(ss);
    ss = ss.(varName{1});
    bouts = ss.bouts;
    btimes = bouts.times(sstates);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calc psd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~, faxis, psd_bouts] = calc_psd('sig', sig, 'bins', btimes,...
    'fs', fs, 'ftarget', ftarget, 'graphics', false);

% interpolate around PLI
if flgPli
    for istate = 1 : nstates
        psd_bouts{istate} = interp_pli(psd_bouts{istate}, faxis);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% remove outlier bouts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% save original before playing
btimesOrig = btimes;
psd_orig = psd_bouts;

% find outliers
clear otl
for istate = 1 : nstates
    data = psd_bouts{istate};
    otl(istate) = get_otl(data, 'thrFactor', 2.5, 'graphics', false);
end

% separate outliers from btimes and psd
otltimes = [];
for istate = 1 : nstates

    idx_bad = otl(istate).idx;

    % times
    otltimes = [otltimes; btimes{istate}(idx_bad, :)];
    btimes{istate}(idx_bad, :) = [];

    % psd
    psd_otl{istate} = psd_bouts{istate}(idx_bad, :);
    psd_bouts{istate}(idx_bad, :) = [];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% organize psd struct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% bouts
psd.bouts.psd = psd_bouts';
psd.bouts.psd_otl = psd_otl;
psd.bouts.times = btimes;
psd.bouts.times_otl = otltimes;

% info
psd.bouts.info.otl = otl;
psd.info.sstates = sstates;
psd.info.snames = snames;
psd.info.clr = clr;
psd.info.faxis = faxis;
psd.info.ch = ch;
psd.info.runtime = datetime("now");
psd.info.flgEmg = flgEmg;
psd.info.flgPli = flgPli;

% calculate mean psd and power in specific frequency bands
for istate = 1 : nstates
    if ~isempty(psd.bouts.psd{istate})

        psd_tmp = psd.bouts.psd{istate};

        % mean
        psd.psd(istate, :) = mean(psd_tmp, 1);

        % band power per bout
        [psd.bands.bouts{istate}, psd.bands.info] = calc_bands('psdData',...
            psd_tmp, 'freq', faxis, 'flgNormBand', false);

        % band power across bouts
        [psd.bands.mean(istate, :), ~] = calc_bands('psdData',...
            psd.psd(istate, :), 'freq', faxis, 'flgNormBand', false);

    end
end

if saveVar
    save(psdfile, 'psd')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~graphics
    return
end

% load data for plotting
if ~exist('emg', 'var')
    emg = load(sleepfile, 'emg_rms');
    emg = emg.emg_rms;
end
load(fullfile(basepath, [basename, '.spec.mat']))

% open figure
setMatlabGraphics(true)
fh = figure;
set(fh, 'WindowState', 'maximized');
tlayout = [3, nstates + 2];
th = tiledlayout(tlayout(1), tlayout(2));
th.TileSpacing = 'tight';
th.Padding = 'none';
title(th, basename, 'interpreter', 'none', 'FontSize', 20)
set(fh, 'DefaultAxesFontSize', 16);

% spectrogram
axh1 = nexttile(th, 1, [1, tlayout(2)]); cla; hold on
plot_spec(spec, 'ch', 1, 'logfreq', false, 'saveFig', false,...
    'axh', axh1, 'xtime', 1)
axis tight
yLimit = ylim;
xticks([]);
xlabel('')

% hypnogram and emg
axh2 = nexttile(th, tlayout(2) + 1, [1, tlayout(2)]); cla; hold on
plot([1 : length(emg)], emg, 'k', 'LineWidth', 0.5)
plot_hypnogram('boutTimes', psd.bouts.times,...
    'clr', clr, 'axh', axh2, 'sstates', [1 : nstates],...
    'yshift', 1)
plot_hypnogram('boutTimes', {psd.bouts.times_otl},...
    'clr', repmat({[0 0 0]}, 1, 1),...
    'axh', axh2, 'sstates', 1, 'yshift', 1.05)
yLimit = ylim;
ylim([yLimit(1), yLimit(2) * 1.05])
xval = [3600 : 3600 : length(emg)];
xticks(xval);
xticklabels(string(xval / 3600))
xlabel('Time [Hr]')
set(axh2, 'YTickMode', 'auto')
set(axh2, 'YColor', 'k')
yticks(axh2, [])
linkaxes([axh1, axh2], 'x')
axis tight
ylabel(axh2, 'EMG')

% psd per state
tilebias = tlayout(2) * 2;
for istate = 1 : nstates
    axh = nexttile(th, istate + tilebias, [1, 1]); cla; hold on

    psdMat = psd.bouts.psd{istate};
    badMat = psd.bouts.psd_otl{istate};
    if ~isempty(badMat)
        ph = plot(faxis, badMat, 'LineWidth', 0.5,...
            'Color', [0.7 0.7 0.7]);
        ph = ph(1);
    end
    if ~isempty(psdMat)
        ph(2) = plot(faxis, mean(psdMat, 1), 'LineWidth', 3,...
            'Color', clr{istate});
    end

    set(gca, 'YScale', 'log', 'XScale', 'log')
    xlabel('Frequency [Hz]')
    ylabel('PSD [mV^2/Hz]')
    title(axh, snames{istate})
    lgdTxt = sprintf('Outliers n=%d', size(badMat, 1));
    legend(ph, {lgdTxt, 'Mean'}, 'Location', 'southwest')
end

% silhouette plot original
tilebias = tlayout(2) * 2 + nstates + 1;
axh = nexttile(th, tilebias, [1, 1]); cla; hold on
plot_silhouette(psd_orig, clr, snames, axh)
title(axh, 'Separation Original')

% silhouette plot after outlier removal
tilebias = tlayout(2) * 2 + nstates + 2;
axh = nexttile(th, tilebias, [1, 1]); cla; hold on
plot_silhouette(psd.bouts.psd, clr, snames, axh)
title(axh, 'Separation Clean')

if flgSaveFig
    figpath = fullfile(basepath, 'graphics', 'sleepState');
    mkdir(figpath)
    figname = fullfile(figpath, [basename, '_', namePrefix]);
    savefig(fh, figname)
end

end

% EOF

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% local functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% interpolate psd to handle power line interference (PLI)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function psd_interp = interp_pli(psd_mat, freqs)

% PLI interpolation preserves ability to analyze full spectrum including gamma
% (30-100 Hz) while removing the confounding 50 Hz line noise artifact. This
% avoids edge effects from split frequency ranges and is more elegant than
% excluding peaks post-hoc. Particularly important when analyzing broadband
% phenomena like cross-frequency coupling. approach adapted from brainstorm.

freq_pli = [48 52];     % frequency range containing line noise [Hz]

% find indices for PLI
idx_pli = freqs >= freq_pli(1) & freqs <= freq_pli(2);
idx_clean = ~idx_pli;

% initialize interpolated psd
psd_interp = psd_mat;

% interpolate across PLI for each spectrum separately
for ispec = 1 : size(psd_mat, 1)
    % get values on either side of PLI
    x = freqs(idx_clean);
    y = psd_mat(ispec, idx_clean);

    % interpolate PLI region
    psd_interp(ispec, idx_pli) = interp1(x, y, freqs(idx_pli), 'pchip');
end

end
