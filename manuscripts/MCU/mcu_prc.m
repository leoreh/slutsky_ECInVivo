clear params
params.winLim = [0 60 * 60];        % Analysis window [s]
params.binSize = 0.001;             % 1ms bins
params.gkHw = 0.012;                % 12ms sigma
params.winStpr = 1.0;               % 1s window
params.nShuffles = 100;            % Number of shuffles
params.spkLim = 2000;

basepaths = [mcu_basepaths('wt_bsl'), mcu_basepaths('mcu_bsl')];
% basepaths = [mcu_basepaths('mea_bac')];
nFiles = length(basepaths);
vars = {'spikes', 'st_metrics', 'fr', 'units'};
v = basepaths2vars('basepaths', basepaths, 'vars', vars);
 

for iFile = 1 : nFiles

    basepath = basepaths{iFile};
    cd(basepath);
    tic
    [prc] = prc_calc(v(iFile).spikes.times, params, 'flgSave', true);
    toc
    prc_plot(prc, 'basepath', basepath, 'flgSave', true);

end




