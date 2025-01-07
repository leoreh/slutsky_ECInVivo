function mea_monoSyn_wrapper(varargin)

% calculates the spike transimission gain between pairs of units.
% Essentially a wrapper for CCH-deconvolution of Lidor.
%
% INPUT:
%   spktimes    cell of spktimes [s].  
%   basepath    char. path to recording session.
%   winCalc     2-element vector describing the window for calculating the
%               stgs [s] {[0 Inf]}.
%   saveVar     logical. save variable {true}.
%   graphics    logical. plot graphics or not {true}. 
%   forceA      logical. reanalyze even if struct file exists
%
% DEPENDENCIES
%   plot_monoSyn
%   CCH-deconvolution (lidor)
%
% TO DO lIST
%   total exc / inh per cell
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
addOptional(p, 'saveVar', true, @islogical);
addOptional(p, 'graphics', true, @islogical);
addOptional(p, 'forceA', false, @islogical);

parse(p, varargin{:})
spktimes    = p.Results.spktimes;
basepath    = p.Results.basepath;
winCalc     = p.Results.winCalc;
saveVar     = p.Results.saveVar;
graphics    = p.Results.graphics;
forceA      = p.Results.forceA;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% params
nunits = length(spktimes);
npairs = nunits * nunits - nunits;
nspks = cellfun(@length, mea.spktimes, 'uni', true);
spkFS = mea.info.fs;

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
for iunit = 1 : nunits
    spkIdx = spktimes{iunit} > winCalc(1) &...
        spktimes{iunit} < winCalc(2);
    spkT = [spkT; spktimes{iunit}(spkIdx)];
    spkL = [spkL; ones(sum(spkIdx), 1) * iunit];
end

% get spk transmission gain. note lidor probably uses a different version
% of CCG.m which requires spktimes in samples.
[excStg, inhStg, act, sil, dcCCH, crCCH, cchbins] =...
    call_cch_stg(spkT * spkFS, spkL, spkFS);

% calc ccg 50 @ 1 counts
[ccg50, ccg50_tstamps] = CCG(mea.spktimes, [], 'binSize', 0.001,...
    'duration',  0.05, 'Fs', 1 / mea.info.fs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% limit results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% limit to pairs with enough counts in the raw cch
ccBad = squeeze(sum(ccg50, 1)) > ccThr;
clear ccExc; [ccExc(1, :), ccExc(2, :)] = find(ccBad);
[~, ccIdx] = cidx2cpair(nunits, [], ccExc');

% limit to pairs with no refractory period
refTidx = round(length(ccg50_tstamps) / 2) - 1 :...
   round(length(ccg50_tstamps) / 2) + 1;
refBad = squeeze(sum(ccg50(refTidx, :, :))) > refThr;
clear refExc; [refExc(1, :), refExc(2, :)] = find(refBad);
[~, refIdx] = cidx2cpair(nunits, [], refExc');

% limit according to stg
stgExcIdx = find(eSTG1 > stgThr(1) | isnan(eSTG1));
stgInhIdx = find(eSTG2 < -stgThr(2) | isnan(eSTG2));

% get remaining excitatory pairs
excIdx = unique([ccIdx, refIdx, stgExcIdx, find(act)]);
excPair = cidx2cpair(nunits, excIdx, []);
fprintf('\nExcitatory: orig = %d, final = %d\n',...
    sum(act), length(excIdx))

% get remaining inhibitory pairs
inhIdx = unique([ccIdx, refIdx, stgInhIdx, find(sil)]);
inhPair = cidx2cpair(nunits, inhIdx, []);
fprintf('\nInhibitory: orig = %d, final = %d\n',...
    sum(sil), length(inhIdx))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if graphics
    
    % excitatory synapses
    npairsPlot = min([nfigs, length(excIdx)]);
    for ipair = 1 : npairsPlot
        cpair = excPair(ipair, :);
        cidx = excIdx(ipair);
        plot_monoSyn('spktimes', mea.spktimes(cpair), 'wv', mea.wv(cpair, :),...
            'wv_std', mea.wv_std(cpair, :), 'clr', 'kk', 'fs', 10000,...
            'saveFig', false, 'units', cpair)
    end
    
    % inhibitory synapses
    npairsPlot = min([nfigs, length(inhIdx)]);
    for ipair = 1 : npairsPlot
        cpair = inhPair(ipair, :);
        cidx = inhIdx(ipair);
        plot_monoSyn('spktimes', mea.spktimes(cpair), 'wv', mea.wv(cpair, :),...
            'wv_std', mea.wv_std(cpair, :), 'clr', 'kk', 'fs', 10000,...
            'saveFig', false, 'units', cpair)
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
monosyn.excPair = excPair;
monosyn.excIdx = excIdx;
monosyn.excStg = excStg;
monosyn.inhPair = inhPair;
monosyn.inhIdx = inhIdx;
monosyn.inhStg = inhStg;

if saveVar
    monofile = fullfile(basepath, [basename, '.monosyn.mat']);
    save(monofile, 'monosyn')
end

end

% EOF

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% legacy messaround
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% 
% % convert a pair of units to an stg idx
% cpair = [14, 17];
% [~, cidx] = cidx2cpair(nunits, [], cpair);
% % convert an stg idx to pair of units
% ccIidx = 11129;
% [cpair, ~] = cidx2cpair(nunits, cidx, []);


% 
% 
% 
% 
% [~, sortidx] = sort(sum(ccg50(:, :, 3), 1));
% cpair = [sortidx(100), 3];
% [~, cidx] = cidx2cpair(nunits, [], cpair);
% plot_monoSyn('spktimes', mea.spktimes(cpair), 'wv', mea.wv(cpair, :),...
%     'wv_std', mea.wv_std(cpair, :), 'clr', 'kk', 'fs', 10000,...
%     'ccg2', dcCCH(:, cidx), 'stg', eSTG2(cidx), 'saveFig', false,...
%     'units', cpair, 'ccg2_tstamps', cchbins * 1000)
% 
% % compute dcCCH. note ccg must be in counts and not rate
% cpair = [1, 2];
% cch = squeeze(ccg50(:, cpair(1), cpair(2)));
% ach11 = squeeze(ccg50(:, cpair(1), cpair(1)));
% ach22 = squeeze(ccg50(:, cpair(2), cpair(2)));
% nspks11 = nspks(cpair(1))
% nspks22 = nspks(cpair(2))
% dcCCH2 = cchdeconv( cch, ach11, nspks11, ach22, nspks22 );
% [~, cidx] = cidx2cpair(nunits, [], cpair);
% max(dcCCH(:, cidx) - dcCCH2)
% 
% [~, pred] = cch_conv( cc( :, u1, u2 ), 1, 'median', 1, 0);
% 
% 
% 
% 
% figure
% subplot(1, 2, 1)
% bar( cchbins, crCCH( :, ccIidx ), 1, 'FaceColor', Cdc, 'EdgeColor', 'none' );
% set( gca, 'box', 'off', 'tickdir', 'out' )
% xlabel( 'Time lag [s]' )
% ylabel( 'Count' )
% title( sprintf( 'rSTG: %0.5f  eSTG: %0.5f', [ NaN eSTG1( ccIidx ) ] ) )
% subplot(1, 2, 2)
% bar( cchbins, dcCCH( :, ccIidx ), 1, 'FaceColor', Cdc, 'EdgeColor', 'none' );
% set( gca, 'box', 'off', 'tickdir', 'out' )
% xlabel( 'Time lag [s]' )
% ylabel( 'Count' )
% title( sprintf( 'rSTG: %0.5f  eSTG: %0.5f', [ NaN eSTG1( ccIidx ) ] ) )
% 
% 

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
