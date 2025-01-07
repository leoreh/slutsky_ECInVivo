function monosyn = monoSyn_wrapper(varargin)

% calculates the spike transimission gain between pairs of units.
% Essentially a wrapper for CCH-deconvolution of Lidor.
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
%   CCH-deconvolution (lidor)
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
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check if already analyzed
[~, basename] = fileparts(basepath);
if ischar(saveVar)
    monofile = fullfile(basepath, [basename, '.monosyn', saveVar, '.mat']);
else
monofile = fullfile(basepath, [basename, '.monosyn.mat']);
end
if exist(monofile) && ~forceA
    load(monofile)
    return
end

% constants
refThr = 1;
ccThr = 500;
stgThr = [0, -Inf];     % [ext, inh]
nfigs = 5;              % number of synapses to plot

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create spktimes and labels in single vec format. limit spks to specified
% window
spkL = [];
spkT = [];
for iunit = 1 : length(spktimes)
    spkIdx = spktimes{iunit} > winCalc(1) &...
        spktimes{iunit} < winCalc(2);
    if sum(spkIdx) == 0
        spkIdx = 1;
    end
    spkT = [spkT; spktimes{iunit}(spkIdx)];
    spkL = [spkL; ones(sum(spkIdx), 1) * iunit];
end
nunits = length(unique(spkL));
npairs = nunits * nunits - nunits;
nspks = cellfun(@length, spktimes, 'uni', true);

% get spk transmission gain. note lidor probably uses a different version
% of CCG.m which requires spktimes in samples.
[eStg, iStg, act, sil, dcCCH, crCCH, cchbins] =...
    call_cch_stg(spkT * fs, spkL, fs);

% calc ccg 50 @ 1 counts
[ccg50, ccg50_tstamps] = CCG(spktimes, [], 'binSize', 0.001,...
    'duration',  0.05, 'Fs', 1 / fs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% limit results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% limit to pairs with enough counts in the raw cch
ccGood = squeeze(sum(ccg50, 1)) > ccThr;
clear ccExc; [ccExc(1, :), ccExc(2, :)] = find(ccGood);
[~, ccIdx] = cidx2cpair(nunits, [], ccExc');

% limit to pairs with no refractory period
refTidx = round(length(ccg50_tstamps) / 2) - 1 :...
   round(length(ccg50_tstamps) / 2) + 1;
refGood = squeeze(sum(ccg50(refTidx, :, :))) > refThr;
clear refExc; [refExc(1, :), refExc(2, :)] = find(refGood);
[~, refIdx] = cidx2cpair(nunits, [], refExc');
cmnIdx = intersect(ccIdx, refIdx);

% limit according to stg
stgEidx = find(eStg > stgThr(1) | ~isnan(eStg));
stgIidx = find(iStg < -stgThr(2) | ~isnan(iStg));

% get remaining excitatory pairs
eIdx = intersect(intersect(cmnIdx, stgEidx), find(act));
ePair = cidx2cpair(nunits, eIdx, []);
fprintf('\nExcitatory: orig = %d, final = %d\n',...
    sum(act), length(eIdx))

% get remaining inhibitory pairs
iIdx = intersect(intersect(cmnIdx, stgIidx), find(sil));
iPair = cidx2cpair(nunits, iIdx, []);
fprintf('\nInhibitory: orig = %d, final = %d\n',...
    sum(sil), length(iIdx))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% synapse stats
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iunit = 1 : nunits
        
    % get indices to all synapses of a unit
    cpair_unit = [ones(1, nunits) * iunit; 1 : nunits]';
    cpair_unit(iunit, :) = []; 
    [~, cidx_unit] = cidx2cpair(nunits, [], cpair_unit);
    
    eNsyn(iunit) = sum(ePair(:, 1) == iunit);
    iNsyn(iunit) = sum(iPair(:, 1) == iunit);
    
    eStrength(iunit) = sum(eStg(intersect(eIdx, cidx_unit)), 'omitnan');
    iStrength(iunit) = sum(iStg(intersect(iIdx, cidx_unit)), 'omitnan');
    
    x = [eNsyn' eStrength'];
    sortrows(x, 2);
       
end

% E-I dominanace (Kobayashi et al., Nat. Comm., 2019)
ei_dom = (eNsyn - iNsyn) ./ (eNsyn + iNsyn);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if graphics
    
    % excitatory synapses
    npairsPlot = min([nfigs, length(eIdx)]);
    for ipair = 1 : npairsPlot
        cpair = ePair(ipair, :);
        cidx = eIdx(ipair);

        plot_monoSyn('spktimes', spktimes(cpair), 'wv', wv(cpair, :),...
            'wv_std', wv_std(cpair, :), 'clr', 'bk', 'fs', 10000,...
            'saveFig', false, 'units', cpair)
    end
    
    % inhibitory synapses
    npairsPlot = min([nfigs, length(iIdx)]);
    for ipair = 1 : npairsPlot
        cpair = iPair(ipair, :);
        cidx = iIdx(ipair);
        plot_monoSyn('spktimes', spktimes(cpair), 'wv', wv(cpair, :),...
            'wv_std', wv_std(cpair, :), 'clr', 'rk', 'fs', 10000,...
            'saveFig', false, 'units', cpair)
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compare manual calculation of dccch to original
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



stms = spktimes(cpair(ipair, :));

% ccg 50 @ 1 counts
[ccg50, ccg50_bins] = CCG(stms, [], 'binSize', 0.001,...
    'duration', 0.05, 'Fs', 1 / fs);
ccg50_bins = ccg50_bins * 1000;


% dcccg 50 @ 1 rate
dcccg = cchdeconv(ccg50(:, 1, 2), ccg50(:, 1, 1), nspks(1),...
    ccg50(:, 2, 2), nspks(2));

% predictor
[~, pred] = cch_conv(dcccg, 11, 'median', 1, 0);

% crcch
dt = diff(ccg50_bins(1 : 2)) / 1000;  % [s]               
den = (ones(length(ccg50_bins), 1) * nspks(1)) * dt; 
crccg = (dcccg - pred) ./ den;

% calculate gain
roiMS = [2, 10]; 
roi_t = ccg50_bins >= roiMS(1) & ccg50_bins <= roiMS(2);
roi_t_idx = find(roi_t);
[g1, g2] = calc_stg(crccg, roi_t, dt, 0, [1, 0]);

% determine if any bin in the ROI is significant
alfa        = 0.001;
nBonf       = sum(roi_t);
gbUpper     = poissinv(1 - alfa / nBonf, max(pred(roi_t, :), [], 1));
gbLower     = poissinv(alfa / nBonf, min( pred( roi_t, : ), [], 1));
act         = any(dcccg(roi_t, :) >= ones(sum(roi_t), 1) * gbUpper, 1);
sil         = any(dcccg(roi_t, :) <= ones(sum(roi_t), 1) * gbLower, 1);
exc_bins    = dcccg(roi_t, :) >= ones(sum(roi_t), 1) * gbUpper;
inh_bins    = dcccg(roi_t, :) <= ones(sum(roi_t), 1) * gbLower;


% ccg 50 @ 1 counts
% [ccg50, ccg50_bins] = CCG(spktimes(cpair), [], 'binSize', 0.001,...
%     'duration', 0.05, 'Fs', 1 / fs);
% ccg50_bins = ccg50_bins * 1000;
% % dcccg 50 @ 1 rate
% dcccg = cchdeconv(ccg50(:, 1, 2), ccg50(:, 1, 1), nspks(cpair(1)),...
%     ccg50(:, 2, 2), nspks(cpair(2)));
% 
% fh = figure;
% subplot(1, 2, 1)
% bh = bar(cchbins, dcCCH(:, cidx), 'BarWidth', 1);
% hold on
% plot([0, 0], ylim, '--k')
% ylabel('Counts')
% xlabel('Time [ms]')
% bh.FaceColor = 'k';
% bh.FaceAlpha = 0.4;
% bh.EdgeColor = 'none';
% box off
% axis tight
% 
% subplot(1, 2, 2)
% bh = bar(cchbins, dcccg, 'BarWidth', 1);
% hold on
% plot([0, 0], ylim, '--k')
% ylabel('Counts')
% xlabel('Time [ms]')
% bh.FaceColor = 'k';
% bh.FaceAlpha = 0.4;
% bh.EdgeColor = 'none';
% box off
% axis tight

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arrange struct and save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

monosyn.info.runtime = datetime(now, 'ConvertFrom', 'datenum');
monosyn.info.basename = basename;
monosyn.info.winCalc = winCalc;
monosyn.info.ccThr = ccThr;
monosyn.info.refThr = refThr;
monosyn.info.stgThr = stgThr;
monosyn.ePair = ePair;
monosyn.eIdx = eIdx;
monosyn.eStg = eStg;
monosyn.iPair = iPair;
monosyn.iIdx = iIdx;
monosyn.iStg = iStg;
monosyn.ei_dom = ei_dom;

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
%     'includeInhibitoryConnections', true, 'bout', [0 4 * 60 * 60]);
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
