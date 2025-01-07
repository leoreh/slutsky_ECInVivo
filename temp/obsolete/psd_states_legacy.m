function psd = psd_states(varargin)

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
%                   time window. can be used, e.g., for baseline period
%   ch              numeric. channels to load from the sigfile. if a vector
%                   is specified, the signal will be averaged across
%                   channels
%   nchans          numeric. no. channels in sigfile
%   fs              numeric. sampling frequency of lfp data. if empty will
%                   be extracted from session struct. also determines the
%                   new sampling frequency of the eeg signal.
%   sigfile         char. name of file to load signal from. if empty but
%                   channel is specified, will load from [basename.lfp]
%   sstates         numeric. index of selected states to calculate psd 
%   ftarget         numeric. requested frequencies for calculating the psd
%   prct            numeric. percent by which to seprate high- and low-emg.
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
% 07 sep 22 LH  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addOptional(p, 'basepath', pwd);
addOptional(p, 'wins', [0 Inf], @isnumeric);
addOptional(p, 'ch', [], @isnumeric);
addOptional(p, 'nchans', [], @isnumeric);
addOptional(p, 'fs', 1, @isnumeric);
addOptional(p, 'sigfile', [], @ischar);
addOptional(p, 'sstates', [1, 4, 5], @isnumeric);
addOptional(p, 'ftarget', [0.5 : 0.5 : 100], @isnumeric);
addOptional(p, 'prct', [70], @isnumeric);
addOptional(p, 'flgEmg', false, @islogical);
addOptional(p, 'saveVar', true, @islogical);
addOptional(p, 'forceA', false, @islogical);
addOptional(p, 'graphics', true, @islogical);

parse(p, varargin{:})
basepath        = p.Results.basepath;
wins            = p.Results.wins;
ch              = p.Results.ch;
nchans          = p.Results.nchans;
fs              = p.Results.fs;
sigfile         = p.Results.sigfile;
sstates         = p.Results.sstates;
ftarget         = p.Results.ftarget;
prct            = p.Results.prct;
flgEmg          = p.Results.flgEmg;
saveVar         = p.Results.saveVar;
forceA          = p.Results.forceA;
graphics        = p.Results.graphics;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% files
cd(basepath)
[~, basename] = fileparts(basepath);
sleepfile = fullfile(basepath, [basename, '.sleep_sig.mat']);

% state params
cfg = as_loadConfig();
if isempty(sstates)
    sstates = 1 : nstates;
end

% load session vars
varsFile = ["sleep_states"; "datInfo"; "session"];
varsName = ["ss"; "datInfo"; "session"];
v = getSessionVars('basepaths', {basepath}, 'varsFile', varsFile,...
    'varsName', varsName);

% check if states file exists, if not flag emg
if isempty(v.ss) || flgEmg
    flgEmg = true;
    psdfile = fullfile(basepath, [basename, '.psdEmg.mat']);
    sstates = [1, 4];    
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

% state bout duration limits
minDur = 10;
interDur = 4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load signal from sSig or from binary if ch specified
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    
    sigfile = fullfile(basepath, [basename, '.sleep_sig.mat']);
    sig = load(sigfile, 'eeg');
    sig = sig.eeg;
    load(sigfile, 'fs');
end


%%% GET ARTIFACTS FROM SPECTROGRAM AND REMOVE FROM EMG IDX


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calc psd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize
% psd.psd = nan(nwin, length(sstates), length(faxis));
for iwin = 1 : nwin
    
    idxWin = floor(wins(iwin, 1) : wins(iwin, 2));

    if ~flgEmg
        boutTimes = v.ss.boutTimes;

    else
        % get indices to high- and low-emg
        labels = double(emg > prctile(emg(idxWin), prct));
        labels(emg < prctile(emg(idxWin), 100 - prct)) = 2;

        % limit indices to time window and get "state" bouts
        labels = labels(idxWin);
        bouts = as_bouts('labels', labels,...
            'minDur', 10, 'interDur', 4);
        boutTimes = bouts.times([1 : 2]);

    end

    % get indices to signal according to window
    idxSig = floor(wins(iwin, :)) * fs;
    if wins(iwin, 1) == 1
        idxSig(1) = 1;
    end
    if idxSig(2) > length(sig)
        idxSig(2) = length(sig);
    end

    % calc psd
    [psd.psd(iwin, :, :), faxis, psd.psd_bouts] = calc_psd('sig',...
        sig(idxSig(1) : idxSig(2), :), 'bins', boutTimes,...
        'fs', fs, 'ftarget', ftarget, 'graphics', graphics);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% organize in struct and save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
psd.info.runtime = datetime(now, 'ConvertFrom', 'datenum');
psd.info.input = p.Results;
psd.info.sigfile = sigfile;
psd.info.fs = fs;
psd.info.wins = wins;
psd.info.faxis = faxis;
psd.info.ftarget = ftarget;
psd.info.sstates = sstates;

if saveVar
    save(psdfile, 'psd')
end    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% graphics flags
flgSmooth = false;
flgNorm = false;
flgSaveFig = true;

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