function psdBins = psd_timebins(varargin)

% calculating the psd per state in different time bins of a
% session.
%
% INPUT:
%   basepath        string. path to recording folder {pwd}
%   fs              numeric. sampling frequency of lfp data. if empty will
%                   be extracted from session struct. also determines the
%                   new sampling frequency of the eeg signal.
%   fsEeg           numeric. sampling frequency of eeg data
%   chEeg           numeric. channels of eeg signal. if greater than 2
%                   will assume it is recorded in same lfp file. otherwise
%                   will assume it is recorded in emg.dat file. if empty
%                   will assume no eeg recording exists. 
%   nchansEeg       numeric. number of channels in [basename].emg.dat file
%   timebins        numeric. timebins to calc psd
%                   in relation to these poitns [sec].  
%   saveVar         logical. save ss var {true}
%   forceA          logical. reanalyze recordings even if ss struct
%                   exists (false)
%   graphics        logical. plot confusion chart and state separation {true}
%
% DEPENDENCIES
%
% TO DO LIST
%
% 07 jan 22 LH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addOptional(p, 'basepath', pwd);
addOptional(p, 'fsEeg', [], @isnumeric);
addOptional(p, 'chEeg', 1, @isnumeric);
addOptional(p, 'nchansEeg', [], @isnumeric);
addOptional(p, 'timebins', [0 Inf], @isnumeric);
addOptional(p, 'tbins_txt', []);
addOptional(p, 'saveVar', true, @islogical);
addOptional(p, 'forceA', false, @islogical);
addOptional(p, 'graphics', true, @islogical);

parse(p, varargin{:})
basepath        = p.Results.basepath;
fsEeg           = p.Results.fsEeg;
chEeg           = p.Results.chEeg;
nchansEeg       = p.Results.nchansEeg;
timebins        = p.Results.timebins;
tbins_txt       = p.Results.tbins_txt;
saveVar         = p.Results.saveVar;
forceA          = p.Results.forceA;
graphics        = p.Results.graphics;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% file
cd(basepath)
[~, basename] = fileparts(basepath);
psdfile = fullfile(basepath, [basename, '.psd_bins.mat']);

% load session vars
varsFile = ["sleep_states"; "datInfo"; "session"];
varsName = ["ss"; "datInfo"; "session"];
v = getSessionVars('basepaths', {basepath}, 'varsFile', varsFile,...
    'varsName', varsName);

% recording params
fileinfo = dir([basename, '.dat']);
nchans = v.session.extracellular.nChannels;
fsDat = v.session.extracellular.sr;
fsLfp = v.session.extracellular.srLfp;
if isempty(fsEeg)
    fsEeg = v.ss.info.sSig.emgFs;
end
if isempty(nchansEeg)
    nchansEeg = nchans;
end
chLfp = v.ss.info.sSig.eegCh;
recLen = floor(fileinfo.bytes / 2 / nchans / fsDat);

nbins = size(timebins, 1);

if isempty(tbins_txt) && nbins == 4
    tbins_txt = {'0-6 ZT', '6-12 ZT', '12-18 ZT', '18-24 ZT'};
else
    tbins_txt = split(num2str(1 : nbins));
end

% params for filtering eeg
import iosr.dsp.*
filtRatio = 450 / (fsEeg / 2);
fsRatio = (fsEeg / fsLfp);

% state params
faxis = 0.2 : 0.2 : 120;
sstates = [1, 4, 5];

% initialize
psdLfp = nan(nbins, length(faxis));
psdEeg = nan(nbins, length(faxis));

% fft params
faxis = [0.2 : 0.2 : 120];
win = hann(2 ^ (nextpow2(2 * fsLfp) - 1));
noverlap = floor(0.25 * fsLfp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calc psd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check if already analyzed
if exist(psdfile, 'file') && ~forceA
    load(psdfile)
else
    for iwin = 1 : nbins
                
        % ---------------------------------------------------------------------
        % hippocampal lfp
        sig = double(bz_LoadBinary([basename, '.lfp'],...
            'duration', diff(timebins(iwin, :)) + 1,...
            'frequency', fsLfp, 'nchannels', nchans, 'start', timebins(iwin, 1),...
            'channels', chLfp, 'downsample', 1));
        sig = mean(sig, 2);
        
        psdLfp(iwin, :) = pwelch(sig, win, noverlap, faxis, fsLfp);
        
        % ---------------------------------------------------------------------
        % frontal eeg
        if ~isempty(chEeg)
            if chEeg < 3
                sig = double(bz_LoadBinary([basename, '.emg.dat'],...
                    'duration', diff(timebins(iwin, :)) + 1,...
                    'frequency', fsEeg, 'nchannels', 2, 'start', timebins(iwin, 1),...
                    'channels', chEeg, 'downsample', 1));
                sig = [iosr.dsp.sincFilter(sig, filtRatio)]';
                sig = sig(fsRatio : fsRatio : length(sig));
            else
                sig = double(bz_LoadBinary([basename, '.lfp'],...
                    'duration', diff(timebins(iwin, :)) + 1,...
                    'frequency', fsLfp, 'nchannels', nchansEeg, 'start', timebins(iwin, 1),...
                    'channels', chEeg, 'downsample', 1));
            end
            psdEeg(iwin, :) = pwelch(sig, win, noverlap, faxis, fsLfp);

        end
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % arrange in struct and save
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
  
    
    psdBins.info.runtime = datetime(now, 'ConvertFrom', 'datenum');
    psdBins.info.chLfp = chLfp;
    psdBins.info.chEeg = chEeg;
    psdBins.info.fsLfp = fsLfp;
    psdBins.info.fsEeg = fsEeg;
    psdBins.timebins = timebins;
    psdBins.psdLfp = psdLfp;
    psdBins.psdEeg = psdEeg;
    
    if saveVar
        save(psdfile, 'psdBins')
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if graphics
    
    if all(isnan(psdBins.psdEeg), 'all')
        nsub = 1;
    else
        nsub = 2;
    end

    setMatlabGraphics(false)
    alphaIdx = linspace(0.5, 1, nbins);
    lim_fAxis = faxis > 1;
    smf = 7;
    gk = gausswin(smf);
    gk = gk / sum(gk);
    
    fh = figure;
    for istate = 1 : length(sstates)
        
        % lfp
        subplot(1, nsub, 1)
        psdMat = squeeze(psdBins.psdLfp(:, lim_fAxis));
        %     for iwin = 1 : nbins
        %         psdMat(iwin, :) = conv(psdMat(iwin, :), gk, 'same');
        %     end
        psdMat = psdMat ./ sum(psdMat, 2);
        ph = plot(faxis(lim_fAxis), psdMat', 'LineWidth', 2);
        for iwin = 1 : length(ph)
            ph(iwin).Color(2) = 1 - iwin * 0.2;
            ph(iwin).Color(3) = 0.5 + iwin * 0.05;
            ph(iwin).Color(4) = alphaIdx(iwin);
        end
        xlabel('Frequency [Hz]')
        ylabel('LFP PSD [mV^2/Hz]')
        legend(tbins_txt, 'Location', 'Southwest')
        set(gca, 'YScale', 'log', 'XScale', 'log')
        xlim([faxis(find(lim_fAxis, 1)), faxis(find(lim_fAxis, 1, 'last'))])
        
        % eeg
        if ~isempty(chEeg)
            subplot(1, nsub, 2)
            psdMat = squeeze(psdBins.psdEeg(:, lim_fAxis));
            %     for iwin = 1 : nbins
            %         psdMat(iwin, :) = conv(psdMat(iwin, :), gk, 'same');
            %     end
            psdMat = psdMat ./ sum(psdMat, 2);
            ph = plot(faxis(lim_fAxis), psdMat', 'LineWidth', 2);
            for iwin = 1 : length(ph)
                ph(iwin).Color(2) = 1 - iwin * 0.2;
                ph(iwin).Color(3) = 0.5 + iwin * 0.05;
                ph(iwin).Color(4) = alphaIdx(iwin);
            end
            xlabel('Frequency [Hz]')
            ylabel('EEG PSD [mV^2/Hz]')
            set(gca, 'YScale', 'log', 'XScale', 'log')
            xlim([faxis(find(lim_fAxis, 1)), faxis(find(lim_fAxis, 1, 'last'))])
        end
    end
    sgtitle(basename)
    
    saveFig = true;
    if saveFig
        figpath = fullfile(basepath, 'graphics');
        mkdir(figpath)
        figname = fullfile(figpath, [basename, '_psdBins']);
        export_fig(figname, '-jpg', '-transparent', '-r300')
    end    
end

