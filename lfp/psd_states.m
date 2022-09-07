function psd = psd_states(varargin)

% wrapper for calc_psd dedicated to sleep states. uses the eeg signal from
% sSig by default, but can load any specified ch from a binary file. can
% separate the recording to windows and calc the psd per state per window.
%
% INPUT:
%   basepath        string. path to recording folder {pwd}
%   wins            numeric. time windows for calculating the psd [sec].
%                   the psd will be calulcated for each state in every
%                   time window
%   ch              numeric. channels to load from the sigfile. if a vector
%                   is specified, the signal will be averaged across
%                   channels
%   fs              numeric. sampling frequency of lfp data. if empty will
%                   be extracted from session struct. also determines the
%                   new sampling frequency of the eeg signal.
%   sigfile         char. name of file to load signal from. if empty but
%                   channel is specified, will load from [basename.lfp]
%   sstates         numeric. index of selected states to calculate psd 
%   faxis           numeric. frequencies for calculating the psd
%   saveVar         logical. save ss var {true}
%   forceA          logical. reanalyze recordings even if ss struct
%                   exists (false)
%   graphics        logical. plot confusion chart and state separation {true}
%
% DEPENDENCIES:
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
addOptional(p, 'fs', 1, @isnumeric);
addOptional(p, 'sigfile', [], @ischar);
addOptional(p, 'sstates', [1, 4, 5], @isnumeric);
addOptional(p, 'faxis', [0.5 : 0.5 : 120], @isnumeric);
addOptional(p, 'saveVar', true, @islogical);
addOptional(p, 'forceA', false, @islogical);
addOptional(p, 'graphics', true, @islogical);

parse(p, varargin{:})
basepath        = p.Results.basepath;
wins            = p.Results.wins;
ch              = p.Results.ch;
fs              = p.Results.fs;
sigfile         = p.Results.sigfile;
sstates         = p.Results.sstates;
faxis           = p.Results.faxis;
saveVar         = p.Results.saveVar;
forceA          = p.Results.forceA;
graphics        = p.Results.graphics;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% files
cd(basepath)
[~, basename] = fileparts(basepath);
psdfile = fullfile(basepath, [basename, '.psd.mat']);

% check if already analyzed
if exist(psdfile, 'file') && ~forceA
    load(psdfile)
    return
end

% load session vars
varsFile = ["sleep_states"; "datInfo"; "session"];
varsName = ["ss"; "datInfo"; "session"];
v = getSessionVars('basepaths', {basepath}, 'varsFile', varsFile,...
    'varsName', varsName);

% assure wins is in range of the recording. note wins is used as index to state labels
if isempty(wins)
    wins = [0 Inf];
end
nwin = size(wins, 1);
wins(wins == 0) = 1;    
wins(wins > length(v.ss.labels)) = length(v.ss.labels);

% state params
cfg = as_loadConfig();
if isempty(sstates)
    sstates = 1 : nstates;
end

% smoothing params (graphics only)
smf = 17;
gk = gausswin(smf);
gk = gk / sum(gk);

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
    sig = mean(sig, 2);

else                    % from sSig
    
    sigfile = fullfile(basepath, [basename, '.sleep_sig.mat']);
    sig = load(sigfile, 'eeg');
    sig = sig.eeg;
    load(sigfile, 'fs');
end

% filter signal? usfeul for emg.dat files
% params for filtering eeg
% import iosr.dsp.*
% filtRatio = 450 / (fsEeg / 2);
% fsRatio = (fsEeg / fsLfp);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calc psd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize
psd.psd = nan(nwin, length(sstates), length(faxis));
for iwin = 1 : nwin

    % calc state epochs in window
    labels = v.ss.labels(wins(iwin, 1) : wins(iwin, 2));
    [stateEpochs, epStats(iwin)] = as_epochs('labels', labels);

    % get indices to signal according to window
    sigidx = wins(iwin, :) * fs;
    if wins(iwin, 1) == 1
        sigidx(1) = 1;
    end
    if sigidx(2) > length(sig)
        sigidx(2) = length(sig);
    end

    % calc psd
    psd.psd(iwin, :, :) = calc_psd('sig', sig(sigidx(1) : sigidx(2)),...
        'bins', stateEpochs(sstates),...
        'fs', fs, 'ftarget', faxis, 'graphics', false);

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