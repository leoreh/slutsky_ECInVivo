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
% 31 jan 22 LH  updates:
% 10 may 25         separated graphics

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
addOptional(p, 'saveVar', true);
addOptional(p, 'graphics', true, @islogical);
addOptional(p, 'saveFig', true, @islogical);
addOptional(p, 'forceA', false, @islogical);

parse(p, varargin{:})
spktimes    = p.Results.spktimes;
basepath    = p.Results.basepath;
winCalc     = p.Results.winCalc;
fs          = p.Results.fs;
saveVar     = p.Results.saveVar;
graphics    = p.Results.graphics;
saveFig     = p.Results.saveFig;
forceA      = p.Results.forceA;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

roi_ms = [1.5, 4];       % region of interest for mono synaptic connections [ms]
alfa = 0.001;            % significance level (consider differentiating for E and I)
refThr = 1;
ccThr = [0 800];         % minimum counts of spikes in cc [E, I]
stgThr = [0.001, -0];    % minimum synpatic strength [E, I]
minConsecBins = [2, 1];  % minimum no of consecutive significant bins [E, I]
nfigs = 20;              % number of synapses to plot

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
[cc.cc50, cc.cc50bins] = CCG(spktimes, [], 'binSize', 0.0005,...
    'duration',  0.05, 'Fs', 1 / fs);
cc.cc50bins = cc.cc50bins * 1000;

% ccg 150 @ 1 counts
[cc.cc150, cc.cc150bins] = CCG(spktimes, [], 'binSize', 0.001,...
    'duration', 0.15, 'Fs', 1 / fs);
cc.cc150bins = cc.cc150bins * 1000;

% ccg 20 @ 0.2 counts
[cc.cc20, cc.cc20bins] = CCG(spktimes, [], 'binSize', 0.0002,...
    'duration', 0.02, 'Fs', 1 / fs);
cc.cc20bins = cc.cc20bins * 1000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calc stg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% run all pairs through nested loops
roi_idx = cc.cc50bins >= roi_ms(1) & cc.cc50bins <= roi_ms(2);
roi_binIdx = find(roi_idx);
roi_nbins = sum(roi_idx);
dt = diff(cc.cc50bins(1 : 2)) / 1000;  % [s]

% initialize
eStg = nan(nunits, nunits);
iStg = nan(nunits, nunits);
eBins = false(roi_nbins, nunits, nunits);
iBins = false(roi_nbins, nunits, nunits);
dccc = zeros(size(cc.cc50));
pred = zeros(size(cc.cc50));
for u1 = 1 : nunits
    for u2 = 1 : nunits
        
        if u1 == u2
            continue
        end
        
        % dccc
        dccc(:, u1, u2) = cchdeconv(cc.cc50(:, u1, u2), cc.cc50(:, u1, u1), nspks(u1),...
            cc.cc50(:, u2, u2), nspks(u2));
        
        % predictor
        [~, pred(:, u1, u2)] = cch_conv(dccc(:, u1, u2), 11, 'median', 1, 0);
        
        % crcch (normalized to rate)
        den = (ones(length(cc.cc50bins), 1) * nspks(u1)) * dt;
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
ccGoode = squeeze(sum(cc.cc50, 1)) > ccThr(1);
ccGoodi = squeeze(sum(cc.cc50, 1)) > ccThr(2);

% limit to synapses with no refractory period
refIdx = round(length(cc.cc50bins) / 2) - 1 :...
   round(length(cc.cc50bins) / 2) + 1;
refGood = squeeze(sum(cc.cc50(refIdx, :, :))) > refThr;

% limit to synapses with minimum stg
eStgGood = eStg > stgThr(1) & ~isnan(eStg);
iStgGood = iStg < stgThr(2) & ~isnan(iStg);

% limit to synapses with minimum significant bins. can also use strfind
% to limit results to consectutive bins.
eConsecGood = false(nunits, nunits);
iConsecGood = false(nunits, nunits);
for u1 = 1 : nunits
    for u2 = 1 : nunits
        
        if u1 == u2
            continue
        end
        eConsecGood(u1, u2) = any(strfind(eBins(:, u1, u2)', ones(1, minConsecBins(1))));
        iConsecGood(u1, u2) = any(strfind(iBins(:, u1, u2)', ones(1, minConsecBins(2))));
    end
end

% get remaining excitatory / inhibitory synapses
eSig = ccGoode & refGood & eStgGood & eConsecGood & squeeze(any(eBins));
iSig = ccGoodi & refGood & iStgGood & iConsecGood & squeeze(any(iBins));
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
% arrange struct and save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[monosyn.eIdx(:, 1), monosyn.eIdx(:, 2)] = find(monosyn.eSig);
[monosyn.iIdx(:, 1), monosyn.iIdx(:, 2)] = find(monosyn.iSig);
monosyn.eBins = eBins;
monosyn.iBins = iBins;
monosyn.eSig = eSig;
monosyn.iSig = iSig;
monosyn.eStg = eStg;
monosyn.iStg = iStg;
monosyn.eiDom = eiDom;
monosyn.nspks = nspks;
monosyn.cc = cc;
monosyn.cc.dccc = dccc;

monosyn.info.runtime = datetime(now, 'ConvertFrom', 'datenum');
monosyn.info.basename = basename;
monosyn.info.winCalc = winCalc;
monosyn.info.ccThr = ccThr;
monosyn.info.refThr = refThr;
monosyn.info.stgThr = stgThr;
monosyn.info.roi_ms = roi_ms;   
monosyn.info.roi_binIdx = roi_binIdx;     
monosyn.info.alfa = alfa;
monosyn.info.pred = pred;
monosyn.info.minConsecBins = minConsecBins;

if saveVar
    save(monofile, 'monosyn')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if graphics
    plot_stg(monosyn, 'basepath', basepath, 'nSyn', 10,...
        'synIdx', [], 'saveFig', true)
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
