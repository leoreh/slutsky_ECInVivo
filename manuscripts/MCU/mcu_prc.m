

clear params
params.winLim = [0 Inf];        % Analysis window [s]
params.binSize = 0.001;             % 1ms bins
params.winStpr = 1.0;               % 1s window
params.nShuffles = 50;             % Number of shuffles
params.spkLim = Inf;

basepaths = [mcu_basepaths('wt_bsl'), mcu_basepaths('mcu_bsl')];
nFiles = length(basepaths);
vars = {'spikes', 'st_metrics', 'fr', 'units', 'spktimes', 'session', 'sleep_states'};
v = basepaths2vars('basepaths', basepaths, 'vars', vars);

idxState = 4; % NREM
recDur = 30 * 60;

for iFile = 1 : nFiles

    % File
    basepath = basepaths{iFile};
    cd(basepath);

    % State
    bouts = v(iFile).ss.bouts.times{idxState};

    fs = v(iFile).session.extracellular.sr;
    mutimes = cellfun(@(x) x / fs, v(iFile).spktimes, 'uni', false);
    spktimes = v(iFile).spikes.times';

    % Restrict to bouts
    spktimes = spktimes_clip(spktimes, bouts, recDur);
    mutimes = spktimes_clip(mutimes, bouts, recDur);

    tic
    prc = prc_calc(v(iFile).spikes.times, ...
        params, ...
        'mutimes', mutimes, ...
        'flgSave', false);
    toc

    rs = mean(prc.prc0_norm(v(iFile).units.type == 'RS'), 'omitnan')
    fs = mean(prc.prc0_norm(v(iFile).units.type == 'FS'), 'omitnan')

    rs = median(prc.pk_norm(v(iFile).units.type == 'RS'), 'omitnan')
    fs = median(prc.pk_norm(v(iFile).units.type == 'FS'), 'omitnan')

    prc_plot(prc, 'basepath', basepath, 'flgSave', false);

end




%% ========================================================================
%  MEA
%  ========================================================================

clear params
recDur = 70 * 60;
params.winLim = [0 recDur];         % Analysis window [s]
params.binSize = 0.001;             % 1ms bins
params.winStpr = 1.0;               % 1s window
params.nShuffles = 50;              % Number of shuffles
params.spkLim = Inf;

basepaths = [mcu_basepaths('mea_bac'), mcu_basepaths('mea_mcuko')];

nFiles = length(basepaths);
vars = {'mea'};
v = basepaths2vars('basepaths', basepaths, 'vars', vars);


for iFile = 2 : nFiles

    % File
    basepath = basepaths{iFile};
    cd(basepath);

    tic
    prc = prc_calc(v(iFile).mea.spktimes, ...
        params, ...
        'flgSave', true);
    toc

    prc_plot(prc, 'basepath', basepath, 'flgSave', true);

end

