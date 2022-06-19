function monosyn = monoSyn_wrapper(varargin)

% calculates the spike transimission gain between pairs of units.
% essentially a wrapper for CCH-deconvolution of Lidor. main differences
% are detection params, graphics, stats, criteria for excluding CCs,
% and here everything runs in a loop s.t. outputs are arranged in a 3d mat
% similar to ccg.m. also ccg input is cell array and not concat vecs
%
% INPUT:
%   spktimes    cell of spktimes [s].  
%   basepath    char. path to recording session.
%   winCalc     2 x n vec describing the window/s for calculating the
%               stgs [s] {[0 Inf]}.
%   fs          numeric. sampling frequency {10000}.
%   wv          2 x n numeric. mean waveform of units.
%   wv_std      2 x n numeric. std of units waveforms.
%   saveVar     logical / char. save variable {true}. if char then save
%               name will include saveVar as a suffix
%   graphics    logical. plot graphics or not {true}. 
%   saveFig     logical. save figure {true}
%   forceA      logical. reanalyze even if struct file exists
%
% DEPENDENCIES
%   plot_monoSyn    
%   cchdeconv       lidor
%   cch_conv        lidor
%   calc_stg        lidor
%   CCG             zugaro (buzcode format)
%
% TO DO lIST
%   separate single synapse analysis
%   separate session analysis from graphics
%   graphics 2 gui
%   implement Kubayashi 2019, https://github.com/NII-Kobayashi/GLMCC
%
% 31 jan 22 LH  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% questions for lidor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% (1) what creiteria do you use to exclude CCHs? perhaps should be separate
% for E / I synapses
% (2) does it make sense that when converting to pairs, the output is sorted
% according to the target cell?
% (3) in cch_stg had to change roiMS to [2 10]. the causility limit of 0.5
% does not exclude many things that look like commen input
% (4) our version of CCG.m receives spktimes in [s]
% (5) no function alines
% (6) in the demo you use a "halfsize" of 100 ms but the default is 50 ms
% (7) the default W in cch_conv is 5 but in cch_stg is 11
% (8) could you please explain the difference between global and local
% maximum in calc_stg? regardless, stgMode 0 does not work unless t_roi is
% changed back to logical
% (9) why not calculate cch at higher resolution (0.2 ms bins) and require
% that at least 2 consecutive bins are significant?
% (10) how do you compare an stg in two time points? simple subtraction or
% index? should i only take those synpases that were significant during
% baseline? what can i do with the rest of the information?
% (11) can you direct towards some knowledge of the mechanism behind stg?
% post- pre-synaptic?

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addOptional(p, 'spktimes', {}, @iscell);
addOptional(p, 'basepath', pwd, @ischar);
addOptional(p, 'winCalc', [0 Inf], @isnumeric);
addOptional(p, 'fs', 10000, @isnumeric);
addOptional(p, 'wv', [], @isnumeric);
addOptional(p, 'wv_std', [], @isnumeric);
addOptional(p, 'saveVar', true);
addOptional(p, 'graphics', true, @islogical);
addOptional(p, 'saveFig', true, @islogical);
addOptional(p, 'forceA', false, @islogical);

parse(p, varargin{:})
spktimes    = p.Results.spktimes;
basepath    = p.Results.basepath;
winCalc     = p.Results.winCalc;
fs          = p.Results.fs;
wv          = p.Results.wv;
wv_std      = p.Results.wv_std;
saveVar     = p.Results.saveVar;
graphics    = p.Results.graphics;
saveFig     = p.Results.saveFig;
forceA      = p.Results.forceA;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

roi_ms = [1, 8];       % region of interest for mono synaptic connections [ms]
alfa = 0.001;          % significance level (consider differentiating for E and I)
refThr = 1;
ccThr = [0 800];
stgThr = [0, -Inf];    % [E, I]
nfigs = 5;             % number of synapses to plot

% for static workspace addition
a = [];
b = [];
c = [];
d = [];
e = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check if already analyzed
[~, basename] = fileparts(basepath);
if ischar(saveVar)
    monofile = fullfile(basepath, [basename, '.monoSyn', saveVar, '.mat']);
else
monofile = fullfile(basepath, [basename, '.monoSyn.mat']);
end
if exist(monofile) && ~forceA
    load(monofile)
    return
end

% limit spktimes to winCalc
funh = @(x) x(InIntervals(x, winCalc)); 
spktimes = cellfun(funh, spktimes, 'uni', false);

% add one spike to empty units. this is important so that
% nunits is independent of winCalc
emptyunits = find(cellfun(@isempty, spktimes, 'uni', true));
for iempty = 1 : length(emptyunits)
    spktimes{emptyunits(iempty)} = 1;
end
nunits = length(spktimes);
npairs = nunits * nunits - nunits;
nspks = cellfun(@length, spktimes, 'uni', true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate CCs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% cc 50 @ 0.5 counts
[cc50, cc50bins] = CCG(spktimes, [], 'binSize', 0.0005,...
    'duration',  0.05, 'Fs', 1 / fs);
cc50bins = cc50bins * 1000;

if graphics
    % ccg 150 @ 1 counts
    [cc150, cc150bins] = CCG(spktimes, [], 'binSize', 0.001,...
        'duration', 0.15, 'Fs', 1 / fs);
    cc150bins = cc150bins * 1000;
    
    % ccg 20 @ 0.2 counts
    [cc20, cc20bins] = CCG(spktimes, [], 'binSize', 0.0002,...
        'duration', 0.02, 'Fs', 1 / fs);
    cc20bins = cc20bins * 1000;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calc stg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% run all pairs through nested loops
roi_idx = cc50bins >= roi_ms(1) & cc50bins <= roi_ms(2);
roi_t = find(roi_idx);
roi_nbins = sum(roi_idx);
dt = diff(cc50bins(1 : 2)) / 1000;  % [s]

% initialize
eStg = nan(nunits, nunits);
iStg = nan(nunits, nunits);
eBins = zeros(roi_nbins, nunits, nunits);
iBins = zeros(roi_nbins, nunits, nunits);
dccc = zeros(size(cc50));
pred = zeros(size(cc50));
for u1 = 1 : nunits
    for u2 = 1 : nunits
        
        if u1 == u2
            continue
        end
        
        % dccc
        dccc(:, u1, u2) = cchdeconv(cc50(:, u1, u2), cc50(:, u1, u1), nspks(u1),...
            cc50(:, u2, u2), nspks(u2));
        
        % predictor
        [~, pred(:, u1, u2)] = cch_conv(dccc(:, u1, u2), 11, 'median', 1, 0);
        
        % crcch (normalized to rate)
        den = (ones(length(cc50bins), 1) * nspks(u1)) * dt;
        crcc = (dccc(:, u1, u2) - pred(:, u1, u2)) ./ den;
        
        % stg
        [eStg(u1, u2), iStg(u1, u2)] = calc_stg(crcc, roi_idx, dt, 0, [1, 0]);
        
        % determine if any bin in the roi is significant
        gbUpper = poissinv(1 - alfa / roi_nbins, max(pred(roi_idx, u1, u2), [], 1));
        gbLower = poissinv(alfa / roi_nbins, min(pred(roi_idx, u1, u2), [], 1));
        eBins(:, u1, u2) = dccc(roi_idx, u1, u2) >= ones(sum(roi_idx), 1) * gbUpper;
        iBins(:, u1, u2) = dccc(roi_idx, u1, u2) <= ones(sum(roi_idx), 1) * gbLower;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% limit results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% limit to synapses with enough cc counts
ccGoode = squeeze(sum(cc50, 1)) > ccThr(1);
ccGoodi = squeeze(sum(cc50, 1)) > ccThr(2);

% limit to synapses with no refractory period
refIdx = round(length(cc50bins) / 2) - 1 :...
   round(length(cc50bins) / 2) + 1;
refGood = squeeze(sum(cc50(refIdx, :, :))) > refThr;

% limit to synapses with minimum stg
eStgGood = eStg > stgThr(1) | ~isnan(eStg);
iStgGood = iStg < -stgThr(2) | ~isnan(iStg);

% limit to synapses with minimum significant bins. can also use strfind
% to limit results to consectutive bins.
% any(strfind(iBins(:, a(isyn))', ones(1, minSigBins)))
% minSigBins = 3;
% [a] = squeeze(sum(eBins, 1) > minSigBins);

% get remaining excitatory / inhibitory synapses
eSig = ccGoode & refGood & eStgGood & squeeze(any(eBins));
iSig = ccGoodi & refGood & iStgGood & squeeze(any(iBins));
fprintf('\nExcitatory: orig = %d, final = %d',...
    sum(squeeze(any(eBins)), 'all'), sum(eSig, 'all'))
fprintf('\nInhibitory: orig = %d, final = %d\n\n',...
    sum(squeeze(any(iBins)), 'all'), sum(iSig, 'all'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% synapse stats
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% E-I dominanace (Kobayashi et al., 2019)
eSyn_unit = sum(eSig, 2);
iSyn_unit = sum(iSig, 2);
eiDom = (eSyn_unit - iSyn_unit) ./ (eSyn_unit + iSyn_unit);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if graphics
    
    % excitatory synapses
    nplot = min([nfigs, sum(eSig, 'all')]);
    if nplot > 0
    [u1, u2] = find(eSig, nplot);
    for isyn = 1 : nplot
        plotSyn(u1(isyn), u2(isyn), 'bk')
    end
    end
    % inhibitory synapses
    nplot = min([nfigs, sum(iSig, 'all')]);
    if nplot > 0
    [u1, u2] = find(iSig, nplot);
    for isyn = 1 : nplot
        plotSyn(u1(isyn), u2(isyn), 'rk')
    end   
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arrange struct and save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

monosyn.info.runtime = datetime(now, 'ConvertFrom', 'datenum');
monosyn.info.basename = basename;
monosyn.info.winCalc = winCalc;
monosyn.info.ccThr = ccThr;
monosyn.info.refThr = refThr;
monosyn.info.stgThr = stgThr;
monosyn.info.roi_ms = roi_ms;     
monosyn.info.alfa = alfa;     
monosyn.eSig = eSig;
monosyn.eStg = eStg;
monosyn.iSig = iSig;
monosyn.iStg = iStg;
monosyn.eiDom = eiDom;
monosyn.nspks = nspks;

if saveVar
    save(monofile, 'monosyn')
end

% nested functions --------------------------------------------------------
    function plotSyn(u1, u2, clr)
        
        % plot mono synaptic connection       
        fh = figure;        
        sb1 = subplot(2, 4, 1);     % ac150 unit1
        sb2 = subplot(2, 4, 2);     % cc50 counts
        sb3 = subplot(2, 4, 3);     % dccc
        sb4 = subplot(2, 4, 4);     % ac150 unit2
        sb5 = subplot(2, 4, 5);     % wv unit1
        sb6 = subplot(2, 4, 6);     % cc150
        sb7 = subplot(2, 4, 7);     % ccg25
        sb8 = subplot(2, 4, 8);     % wv unit2
        
        % acg 1
        set(gcf, 'CurrentAxes', sb1)
        plot_ccg(cc150(:, u1, u1), cc150bins, 'clr', clr(1),...
            'pred', [], 'sigbins1', [], 'sigbins2', [])
        title(sprintf('Unit #%d (Presynaptic)', u1))
        
        % acg 2
        set(gcf, 'CurrentAxes', sb4)
        plot_ccg(cc150(:, u2, u2), cc150bins, 'clr', clr(2),...
            'pred', [], 'sigbins1', [], 'sigbins2', [])
        title(sprintf('Unit #%d (Postsynaptic)', u2))
        
        % cc50 counts
        set(gcf, 'CurrentAxes', sb2)
        plot_ccg(cc50(:, u1, u2), cc50bins, 'clr', 'k',...
            'pred', [], 'sigbins1', roi_t(logical(eBins(:, u1, u2))),...
            'sigbins2', roi_t(logical(iBins(:, u1, u2))))
        
        % dccc
        set(gcf, 'CurrentAxes', sb3)
        plot_ccg(dccc(:, u1, u2), cc50bins, 'clr', 'k',...
            'pred', pred(:, u1, u2), 'sigbins1', roi_t(logical(eBins(:, u1, u2))),...
            'sigbins2', roi_t(logical(iBins(:, u1, u2))))
        if strcmp(clr(1), 'r')
            title(sprintf('iSTG = %.4f', iStg(u1, u2)))
        else
            title(sprintf('eSTG = %.4f', eStg(u1, u2)))
        end
        
        % cc150 counts
        set(gcf, 'CurrentAxes', sb6)
        plot_ccg(cc150(:, u1, u2), cc150bins, 'clr', 'k',...
            'pred', [], 'sigbins1', [], 'sigbins2', [])
        
        % cc20 counts
        set(gcf, 'CurrentAxes', sb7)
        plot_ccg(cc20(:, u1, u2), cc20bins, 'clr', 'k',...
            'pred', [], 'sigbins1', [], 'sigbins2', [])
        
        if ~isempty(wv)
            
            % waveform 1
            set(gcf, 'CurrentAxes', sb5)
            x_val = [1 : size(wv, 2)] / fs * 1000;
            plot(x_val, wv(u1, :), clr(1), 'LineWidth', 2)
            if ~isempty(wv_std)
                patch([x_val, flip(x_val)], [wv(u1, :) + wv_std(u1, :),...
                    flip(wv(u1, :) - wv_std(u1, :))],...
                    clr(1), 'EdgeColor', 'none', 'FaceAlpha', .2, 'HitTest', 'off')
            end
            xlabel('Time [ms]')
            ylabel('Voltage [mV]')
            
            % waveform 2
            set(gcf, 'CurrentAxes', sb8)
            x_val = [1 : size(wv, 2)] / fs * 1000;
            plot(x_val, wv(u2, :), clr(2), 'LineWidth', 2)
            if ~isempty(wv_std)
                patch([x_val, flip(x_val)], [wv(u2, :) + wv_std(u2, :),...
                    flip(wv(u2, :) - wv_std(u2, :))],...
                    clr(2), 'EdgeColor', 'none', 'FaceAlpha', .2, 'HitTest', 'off')
            end
            xlabel('Time [ms]')
            ylabel('Voltage [mV]')
            
        end
        drawnow
        
        % save
        if saveFig
            figpath = fullfile(basepath, 'graphics', 'monoSyn');
            figname = fullfile(figpath, sprintf('monoSyn_%d_%d', u1, u2));
            export_fig(figname, '-jpg', '-transparent', '-r300')
        end
        
    end

end

% EOF

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cell explorer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% compare to mono_res
% sortrows(mono_res.sig_con_inhibitory, 2)
% sortrows(mono_res.sig_con_excitatory, 2)


% no need to limit ccg. just do manual. 
% must figure out time stamps of mea
% 
% spikes.times = mea.spktimes;
% spikes.shankID = mea.ch;
% spikes.cluID = [1 : length(mea.spktimes)];
% 
% mono_res = ce_MonoSynConvClick(spikes,...
%     'includeInhibitoryConnections', true, 'epoch', [0 4 * 60 * 60]);
% 
% basepath = pwd;
% [~, basename] = fileparts(basepath);
% monofile = fullfile(basepath, [basename, '.mono_res2.cellinfo.mat']);
% save(monofile, 'mono_res');
% setMatlabGraphics(true)
% gui_MonoSyn(monofile);
% load(monofile)
% mono_res
% 
% cpair = [93, 54];
% 
% clr = 'rb';
% plot_monoSyn('spktimes', mea.spktimes(cpair), 'wv', mea.wv(cpair, :),...
%     'wv_std', mea.wv_std(cpair, :), 'clr', clr, 'fs', 10000,...
%     'ccg2', [], 'stg', [], 'saveFig', false,...
%     'units', cpair, 'ccg2_tstamps', [])
