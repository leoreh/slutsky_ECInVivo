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
%   winCalc     2-element vector describing the window for calculating the
%               stgs [s] {[0 Inf]}.
%   fs          numeric. sampling frequency {10000}.
%   wv          2 x n numeric. mean waveform of units.
%   wv_std      2 x n numeric. std of units waveforms.
%   saveVar     logical / char. save variable {true}. if char then save
%               name will include saveVar as a suffix
%   graphics    logical. plot graphics or not {true}. 
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
%   calc synapse stats per cell. include number of significant bins
%
% 31 jan 22 LH  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% questions for lidor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% (1) what creiteria do you use to exclude CCHs?
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
validate_win = @(win) assert(isnumeric(win) && length(win) == 2,...
    'time window must be in the format [start end]');

p = inputParser;
addOptional(p, 'spktimes', {}, @iscell);
addOptional(p, 'basepath', pwd, @ischar);
addOptional(p, 'winCalc', [0 Inf], validate_win);
addOptional(p, 'fs', 10000, @isnumeric);
addOptional(p, 'wv', [], @isnumeric);
addOptional(p, 'wv_std', [], @isnumeric);
addOptional(p, 'saveVar', true);
addOptional(p, 'graphics', true, @islogical);
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
forceA      = p.Results.forceA;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

roi_ms = [2, 10];       % region of interest for mono synaptic connections [ms]
alfa = 0.001;           % significance level (consider differentiating for E and I)
refThr = 1;
ccThr = 500;
stgThr = [0, -Inf];     % [E, I]
nfigs = 5;              % number of synapses to plot

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

nunits = length(spktimes);
npairs = nunits * nunits - nunits;
nspks = cellfun(@length, spktimes, 'uni', true);

% cc 50 @ 1 counts
[cc50, cc50bins] = CCG(spktimes, [], 'binSize', 0.001,...
    'duration',  0.05, 'Fs', 1 / fs);
cc50bins = cc50bins * 1000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calc stg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% run all pairs through nested loops
roi_idx = cc50bins >= roi_ms(1) & cc50bins <= roi_ms(2);
roi_nbins = sum(roi_idx);
dt = diff(cc50bins(1 : 2)) / 1000;  % [s]

% initialize
eStg = nan(nunits, nunits);
iStg = nan(nunits, nunits);
eBins = zeros(roi_nbins, nunits, nunits);
iBins = zeros(roi_nbins, nunits, nunits);
for u1 = 1 : nunits
    for u2 = 1 : nunits
        
        if u1 == u2
            continue
        end
        
        % dccc
        dccc = cchdeconv(cc50(:, u1, u2), cc50(:, u1, u1), nspks(u1),...
            cc50(:, u2, u2), nspks(u2));
        
        % predictor
        [~, pred] = cch_conv(dccc, 11, 'median', 1, 0);
        
        % crcch (normalized to rate)
        den = (ones(length(cc50bins), 1) * nspks(u1)) * dt;
        crcc = (dccc - pred) ./ den;
        
        % stg
        [eStg(u1, u2), iStg(u1, u2)] = calc_stg(crcc, roi_idx, dt, 0, [1, 0]);
        
        % determine if any bin in the roi is significant
        gbUpper = poissinv(1 - alfa / roi_nbins, max(pred(roi_idx, :), [], 1));
        gbLower = poissinv(alfa / roi_nbins, min( pred( roi_idx, : ), [], 1));
        eBins(:, u1, u2) = dccc(roi_idx, :) >= ones(sum(roi_idx), 1) * gbUpper;
        iBins(:, u1, u2) = dccc(roi_idx, :) <= ones(sum(roi_idx), 1) * gbLower;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% limit results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% limit to synapses with enough cc counts
ccGood = squeeze(sum(cc50, 1)) > ccThr;

% limit to synapses with no refractory period
refIdx = round(length(cc50bins) / 2) - 1 :...
   round(length(cc50bins) / 2) + 1;
refGood = squeeze(sum(cc50(refIdx, :, :))) > refThr;

% limit to synapses with minimum stg
eStgGood = eStg > stgThr(1) | ~isnan(eStg);
iStgGood = iStg < -stgThr(2) | ~isnan(iStg);

% get remaining excitatory / inhibitory synapses
eSig = ccGood & refGood & eStgGood & squeeze(any(eBins));
iSig = ccGood & refGood & iStgGood & squeeze(any(iBins));
fprintf('\nExcitatory: orig = %d, final = %d',...
    sum(squeeze(any(eBins)), 'all'), sum(squeeze(any(eSig)), 'all'))
fprintf('\nInhibitory: orig = %d, final = %d\n\n',...
    sum(squeeze(any(iBins)), 'all'), sum(squeeze(any(iSig)), 'all'))

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
    [u1, u2] = find(eSig, nplot);
    for isyn = 1 : nplot
        usyn = [u1(isyn), u2(isyn)];
        plot_monoSyn('spktimes', spktimes(usyn), 'wv', wv(usyn, :),...
            'wv_std', wv_std(usyn, :), 'clr', 'bk', 'fs', 10000,...
            'saveFig', false, 'units', usyn)
    end
    
    % inhibitory synapses
    nplot = min([nfigs, sum(iSig, 'all')]);
    [u1, u2] = find(iSig, nplot);
    for isyn = 1 : nplot
        usyn = [u1(isyn), u2(isyn)];
        plot_monoSyn('spktimes', spktimes(usyn), 'wv', wv(usyn, :),...
            'wv_std', wv_std(usyn, :), 'clr', 'bk', 'fs', 10000,...
            'saveFig', false, 'units', usyn)
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
