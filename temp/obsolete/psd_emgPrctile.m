function psd = psd_emgPrctile(varargin)

% wrapper for calc_psd that uses the emg_rms from sSig to calculate the psd
% during high- and low-emg activity. uses the eeg signal from sSig by
% default, but can load any specified ch from a binary file.
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
%   prct            numeric. percent by which to seprate high- and low-emg.
%                   e.g., prct = 50 is the median 
%   ftarget         numeric. requested frequencies for calculating the psd
%   saveVar         logical. save var {true}
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
addOptional(p, 'fs', 1, @isnumeric);
addOptional(p, 'sigfile', [], @ischar);
addOptional(p, 'prct', [70], @isnumeric);
addOptional(p, 'ftarget', [0.5 : 0.5 : 120], @isnumeric);
addOptional(p, 'saveVar', true, @islogical);
addOptional(p, 'forceA', false, @islogical);
addOptional(p, 'graphics', true, @islogical);

parse(p, varargin{:})
basepath        = p.Results.basepath;
wins            = p.Results.wins;
ch              = p.Results.ch;
fs              = p.Results.fs;
sigfile         = p.Results.sigfile;
prct            = p.Results.prct;
ftarget         = p.Results.ftarget;
saveVar         = p.Results.saveVar;
forceA          = p.Results.forceA;
graphics        = p.Results.graphics;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% files
cd(basepath)
[~, basename] = fileparts(basepath);
psdfile = fullfile(basepath, [basename, '.psdEmg.mat']);
sleepfile = fullfile(basepath, [basename, '.sleep_sig.mat']);

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

% load emg
emg = load(sleepfile, 'emg_rms');
emg = emg.emg_rms;

% assure wins is in range of recording. note wins is used as index to state labels
if isempty(wins)
    wins = [0 Inf];
end
nwin = size(wins, 1);
wins(wins == 0) = 1;    
wins(wins > length(emg)) = length(emg);

% state params
cfg = as_loadConfig();
sstates = [1, 4];       % AW and NREM

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
    if length(ch) > 1
        sig = mean(sig, 2);
    end

else                    % from sSig
    
    sigfile = fullfile(basepath, [basename, '.sleep_sig.mat']);
    sig = load(sigfile, 'eeg');
    sig = sig.eeg;
    load(sigfile, 'fs');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calc psd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize
% psd.psd = nan(nwin, length(sstates), length(faxis));
for iwin = 1 : nwin
 
    % get indices to high- and low-emg
    labels = double(emg > prctile(emg, prct));
    labels(emg < prctile(emg, 100 - prct)) = 2;
    
    % limit indices to time window and get "state" bouts
    labels = labels(wins(iwin, 1) : wins(iwin, 2));
    [boutTimes, ~] = as_bouts('labels', labels, 'minDur', 10, 'interDur', 4);
    boutTimes = boutTimes([1 : 2]);

    % get indices to signal according to window
    sigidx = wins(iwin, :) * fs;
    if wins(iwin, 1) == 1
        sigidx(1) = 1;
    end
    if sigidx(2) > length(sig)
        sigidx(2) = length(sig);
    end

    % calc psd
    [psd.psd(iwin, :, :), faxis] = calc_psd('sig',...
        sig(sigidx(1) : sigidx(2), :), 'bins', boutTimes,...
        'fs', fs, 'ftarget', ftarget, 'graphics', false);

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

    for istate = 1 : 2      
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
        ylim([0.1, 10^5])
        
    end
    title(th, basename)
    
    if flgSaveFig
        figpath = fullfile(basepath, 'graphics', 'sleepState');
        mkdir(figpath)
        figname = fullfile(figpath, [basename, '_psdEmg']);
        export_fig(figname, '-jpg', '-transparent', '-r300')
    end  
end

% EOF