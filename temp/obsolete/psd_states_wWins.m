function psd = psd_states(varargin)

% NOTE AFTER MULTIPLE CHANGES, WIN FUNCTIONALITY PROBABLY NOT WORKING
% wrapper for calc_psd dedicated to sleep states. uses the eeg signal from
% sSig by default, but can load any specified ch from a binary file. can
% separate the recording to windows and calc the psd per state per window.
% if sleep_states.mat is not found, will separate the recording to "AW" and
% "NREM" according to high- and low-emg activity. 
%
% INPUT:
%   basepath        string. path to recording folder {pwd}
%   wins            numeric. time windows for calculating the psd [sec].
%                   the psd will be calulcated for each state in every
%                   time window. For example, can be used to limit the
%                   analysis to baseline period or to separate light and
%                   dark phases
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
%   boutTimes       cell of n x 2 mats. if empty will calculate from
%                   ss.labels or from emg_labels
%   sstates         numeric. index of selected states to calculate psd 
%   ftarget         numeric. requested frequencies for calculating the psd
%   emgThr          numeric. percent by which to seprate high- and low-emg.
%                   e.g., prct = 50 is the median 
%   flgEmg          logical. calc psd in high- and low-emg even if states
%                   file exists
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addOptional(p, 'basepath', pwd);
addOptional(p, 'wins', [0 Inf], @isnumeric);
addOptional(p, 'sig', [], @isnumeric);
addOptional(p, 'ch', [], @isnumeric);
addOptional(p, 'nchans', [], @isnumeric);
addOptional(p, 'fs', 1, @isnumeric);
addOptional(p, 'sigfile', [], @ischar);
addOptional(p, 'sstates', [1, 4, 5], @isnumeric);
addOptional(p, 'boutTimes', []);
addOptional(p, 'ftarget', [0.5 : 0.5 : 100], @isnumeric);
addOptional(p, 'emgThr', [50], @isnumeric);
addOptional(p, 'flgEmg', false, @islogical);
addOptional(p, 'saveVar', true, @islogical);
addOptional(p, 'forceA', false, @islogical);
addOptional(p, 'graphics', true, @islogical);

parse(p, varargin{:})
basepath        = p.Results.basepath;
wins            = p.Results.wins;
sig             = p.Results.sig;
ch              = p.Results.ch;
nchans          = p.Results.nchans;
fs              = p.Results.fs;
sigfile         = p.Results.sigfile;
sstates         = p.Results.sstates;
boutTimes       = p.Results.boutTimes;
ftarget         = p.Results.ftarget;
emgThr          = p.Results.emgThr;
flgEmg          = p.Results.flgEmg;
saveVar         = p.Results.saveVar;
forceA          = p.Results.forceA;
graphics        = p.Results.graphics;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% graphics flags
flgSmooth = false;
flgNorm = false;
flgSaveFig = true;

% state params
cfg = as_loadConfig();
if isempty(sstates)
    sstates = 1 : nstates;
end

% state bout duration limits
minDur = 20;
interDur = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% file
cd(basepath)
[~, basename] = fileparts(basepath);
sleepfile = fullfile(basepath, [basename, '.sleep_sig.mat']);

% load session vars
varsFile = ["sleep_states"; "datInfo"; "session"];
varsName = ["ss"; "datInfo"; "session"];
v = getSessionVars('basepaths', {basepath}, 'varsFile', varsFile,...
    'varsName', varsName);

% check if states file exists, if not flag emg
if isempty(v.ss) || flgEmg
    flgEmg = true;
    psdfile = fullfile(basepath, [basename, '.psdEmg.mat']);
    sstates = [1, 2];    
    emg = load(sleepfile, 'emg_rms');
    emg = emg.emg_rms;
else
    flgEmg = false;
    psdfile = fullfile(basepath, [basename, '.psd.mat']);
end

% check if already analyzed
if exist(psdfile, 'file') && ~forceA
    load(psdfile)
    return
end

% assure wins is in range of recording. note wins is used as index to state labels
if isempty(wins)
    wins = [0 Inf];
end
nwin = size(wins, 1);
wins(wins == 0) = 1;    
if ~flgEmg
    wins(wins > length(v.ss.labels)) = length(v.ss.labels);
else
    wins(wins > length(emg)) = length(emg);
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
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate state bouts if not provided
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iwin = 1 : nwin
    
    idxWin = floor(wins(iwin, 1) : wins(iwin, 2));

    % get boutTimes
    if isempty(boutTimes)
        if ~flgEmg           
            
            % get AS labels and limit to time window
            labels = v.ss.labels(idxWin);
        else
            
            % find threshold to separate the bimodal distribution of emg
            if isempty(emgThr)
                [~, cents] = kmeans(emg(:), 2);
                emgThr = mean(cents);
            end
            
            % create EMG labels and limit to time window
            labels = double(emg > emgThr);
            labels(emg < emgThr) = 2;
            labels = labels(idxWin);
        end
        
        % re-calc state bouts from labels
        bouts = as_bouts('labels', labels,...
            'minDur', minDur, 'interDur', interDur, 'rmArtifacts', true);
        boutTimes = bouts.times;
    end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calc psd and bands
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iwin = 1 : nwin
    
    % get indices to signal according to window
    idxSig = floor(wins(iwin, :)) * fs;
    if wins(iwin, 1) == 1
        idxSig(1) = 1;
    end
    if idxSig(2) > length(sig)
        idxSig(2) = length(sig);
    end

    % calc psd
    [psd.psd(iwin, :, :), faxis, psd.psd_bouts(iwin, :)] = calc_psd('sig',...
        sig(idxSig(1) : idxSig(2), :), 'bins', boutTimes(sstates),...
        'fs', fs, 'ftarget', ftarget, 'graphics', graphics);
    
    % calc power in specific frequency bands
    for istate = 1 : length(sstates)
        stateIdx = sstates(istate);
        [psd.bands{iwin, istate}, psd.info.bands] = calc_bands('psdData', psd.psd_bouts{istate},...
            'freq', faxis, 'flgNormBand', false);
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% organize in struct and save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
psd.info.boutTimes = bouts.times;
psd.info.runtime = datetime("now");
psd.info.input = p.Results;
psd.info.sigfile = sigfile;
psd.info.fs = fs;
psd.info.wins = wins;
psd.info.faxis = faxis;
psd.info.ftarget = ftarget;
psd.info.sstates = sstates;
psd.info.emgThr = emgThr;

if saveVar
    save(psdfile, 'psd')
end    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if graphics
             
    setMatlabGraphics(false)
    fh = figure;
    fh.Position = [0.1 0.1 0.8 0.8];
    th = tiledlayout(1, length(sstates), 'TileSpacing', 'Compact');

    % plot params
    alphaIdx = linspace(0.5, 1, nwin);
    lim_fAxis = faxis >= 1;

    for istate = 1 : length(sstates)        
        axh = nexttile;
        
        % grab relevant data
        psdMat = squeeze(psd.psd(:, istate, lim_fAxis));
        if isvector(psdMat)
            psdMat = psdMat';
        end

        if flgSmooth
            smf = 17;
            gk = gausswin(smf);
            gk = gk / sum(gk);

            for iwin = 1 : nwin
                psdMat(iwin, :) = conv(psdMat(iwin, :), gk, 'same');
            end
        end

        if flgNorm
            psdMat = psdMat ./ sum(psdMat, 2);
            ytxt = 'Norm. LFP PSD';
        else
            ytxt = 'PSD [mV^2/Hz]';
        end

        ph = plot(faxis(lim_fAxis), psdMat', 'LineWidth', 2);
        for iwin = 1 : length(ph)
            ph(iwin).Color(istate) = cfg.colors{sstates(istate)}(istate) - iwin * 0.003;
            ph(iwin).Color(4) = alphaIdx(iwin);
        end
        set(gca, 'YScale', 'log', 'XScale', 'log')
        title(axh, cfg.names{sstates(istate)})
        xlabel('Frequency [Hz]')
        ylabel(ytxt)
        legend('Location', 'Southwest', 'FontSize', 9)
        xlim([faxis(find(lim_fAxis, 1)), faxis(find(lim_fAxis, 1, 'last'))])
         
    end
    title(th, basename)
    
    if flgSaveFig
        figpath = fullfile(basepath, 'graphics', 'sleepState');
        mkdir(figpath)
        figname = fullfile(figpath, [basename, '_psdStates']);
        export_fig(figname, '-jpg', '-transparent', '-r300')
    end  
end

% EOF