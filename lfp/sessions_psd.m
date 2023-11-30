function [bands, powdb] = sessions_psd(mname, varargin)

% arranges the psd in states from different sessions of the same mouse.
% also calculates the power in specific bands. can normalize the psd to
% broadband power and can normalize each session to a baseline session/s.
% plots the results.
%
% INPUT:
%   mname           string. mouse name to analyze
%   idxBsl          numeric. indices to baseline sessions for normalization
%   flgNormTime     logical. normalize each session to baseline (1st
%                   session) {true}. only applied to bands
%   flgNormBand     logical. normalize psd to broadband {true}
%   flgEmg          logical. load [basename].psdEmg.mat {true}
%   flgAnalyze      logical. re-analyze psd {false}
%   flgDb           logical. convert psd to dB {true}
%   saveVar         logical. save ss var {true}
%   graphics        logical. plot confusion chart and state separation {true}
%
% OUTPUT:
%   bands           numeric 3d mat of band x state x session
%   powdb           numeric 3d mat of state x freq x session
% 
% DEPENDENCIES:
%   calc_psd        
%
% TO DO LIST:
%
% 20 oct 22 LH  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addOptional(p, 'idxBsl', 1, @isnumeric);
addOptional(p, 'flgNormTime', true, @islogical);
addOptional(p, 'flgNormBand', true, @islogical);
addOptional(p, 'flgEmg', true, @islogical);
addOptional(p, 'flgAnalyze', false, @islogical);
addOptional(p, 'flgDb', true, @islogical);
addOptional(p, 'saveVar', true, @islogical);
addOptional(p, 'graphics', true, @islogical);

parse(p, varargin{:})
idxBsl          = p.Results.idxBsl;
flgNormTime     = p.Results.flgNormTime;
flgNormBand     = p.Results.flgNormBand;
flgEmg          = p.Results.flgEmg;
flgAnalyze      = p.Results.flgAnalyze;
flgDb           = p.Results.flgDb;
saveVar         = p.Results.saveVar;
graphics        = p.Results.graphics;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load or analyze
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load vars from each session
if flgEmg
    varsFile = ["psdEmg"];
else
    varsFile = ["psd"];
end
varsName = ["psd"];
[v, basepaths] = getSessionVars('mname', mname, 'varsFile', varsFile,...
    'varsName', varsName);
nfiles = length(basepaths);

% re-analyze
if flgAnalyze
    
    % grab params from first psd file to make sure all PSDs are analyzed
    % the same. otherwise use defaults.
    if ~isempty(v(1).psd)
        ch = v(1).psd.info.input.ch;
        prct = v(1).psd.info.input.prct;
        sstates = v(1).psd.info.input.sstates;
        ftarget = v(1).psd.info.input.ftarget;
    else
        ch = [];
        prct = 70;
        sstates = [1, 4];
        ftarget = [0.5 : 0.5 : 100];
    end

    % go over files and re-analyze
    for ifile = 1 : nfiles
        
        % get time window 
        if ~isempty(v(ifile).psd)
            wins = v(ifile).psd.info.wins;
        else
            wins = [0, Inf];
        end

        % calc psd
        ch = [13 : 16];
        v(ifile).psd = psd_states('basepath', basepaths{ifile}, 'sstates', sstates,...
            'ch', ch, 'fs', [], 'wins', wins, 'saveVar', true,...
            'graphics', false, 'forceA', true, 'ftarget', ftarget,...
            'prct', prct, 'flgEmg', flgEmg);
    end
end

% cat psd struct
psd = catfields([v.psd], 'catdef', 'addim');
powdb = psd.psd;

% get updated params
cfg = as_loadConfig();
sstates = psd.info.sstates(:, 1);
nstates = length(sstates);
freq = v(1).psd.info.faxis;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bands
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% band params. taken from Boyce et al., Science, 2016
bandNames = ["broad", "swa", "delta", "theta", "alpha", "beta", "gamma"];
bandFreqs = [freq(1), 100; 0.5, 1; 1, 4; 4, 10; 10, 14; 15, 30; 30, 100];

% convert to dB (chronux output is power not magnitude)
if flgDb
    powdb = 10 * log10(powdb);
    ytxt = 'PSD [dB]';
else
    ytxt = 'PSD [mV^2/Hz]';
end

% normalize to broadband
lim_freq = freq >= 0 & freq < Inf;
if flgNormBand
    powdb = powdb ./ sum(powdb(:, lim_freq, :), 2);
    sbands = [2 : length(bandNames)];
    ytxt = 'Norm. PSD';
else
    sbands = [1 : length(bandNames)];
end

% calc power in band. bands is a 3d array of freqBand x state x session.
% psd is a 3d array of state x freq x session.
for iband = 1 : length(bandFreqs)
    bandIdx = InIntervals(freq, bandFreqs(iband, :));
    bands(iband, :, :) = squeeze(sum(powdb(:, bandIdx, :), 2));
end

% normalize to baseline
if flgNormTime
    bsl = mean(bands(:, :, idxBsl), 3, 'omitnan');
    bands = (bands ./ bsl) * 100;
else
    bands = 10 * log10(bands);
end

% organize and save
if saveVar
    mousepath = fileparts(basepaths{1});
    bandsfile = fullfile(mousepath, [mname, '_psdBands.mat']);
    save(bandsfile, 'bands')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if graphics

    fh = figure;
    th = tiledlayout(2, length(sstates), 'TileSpacing', 'Compact');
    alphaIdx = linspace(0.5, 1, nfiles);

    for istate = 1 : nstates
        axh = nexttile;

        ydata = squeeze(powdb(istate, :, :));
        ph = plot(freq, ydata, 'LineWidth', 2);
        for ifile = 1 : nfiles
            ph(ifile).Color(istate) = cfg.colors{sstates(istate)}(istate) - ifile * 0.04;
            ph(ifile).Color(4) = alphaIdx(ifile);
        end
        set(gca, 'XScale', 'log')
        if ~flgDb
            set(gca, 'YScale', 'log')
        end
        title(cfg.names{sstates(istate)})
        xlabel('Frequency [Hz]')
        ylabel(ytxt)
        legend(split(num2str(1 : nfiles)), 'Location', 'Southwest', 'FontSize', 9)

    end

    for istate = 1 : nstates
        axh = nexttile;
        plot(squeeze(bands(sbands, istate, :))', 'LineWidth', 2)
        hold on
        plot(xlim, [100 100], '--k')
        ylabel([ytxt, ' (% BSL)'])
        set(gca, 'YScale', 'log')
    end
    legend(bandNames(sbands))
    title(th, mname)

    % save fig
     figpath = fullfile(mousepath, 'graphics');
     mkdir(figpath)
     figname = fullfile(figpath, [mname, '_psdBands']);
     export_fig(figname, '-jpg', '-transparent', '-r300')

end

% EOF