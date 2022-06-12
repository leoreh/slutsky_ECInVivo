

% params
cfg = as_loadConfig();
nstates = cfg.nstates;
sstates = [1, 4];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% per session spec vs. time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

isession = 2;
[expData, xData] = sessions_catVarTime('mname', '',...
    'dataPreset', {'fr', 'sleep_emg', 'spec', 'bands'}, 'graphics', true,...
    'basepaths', basepaths(isession), 'xTicksBinsize', 3,...
    'markRecTrans', true);

% spec
spec = calc_spec('sig', [], 'fs', 1250, 'graphics', true, 'saveVar', false,...
    'padfft', -1, 'winstep', 5, 'logfreq', false, 'ftarget', [],...
    'ch', v(isession).spec.info.ch, 'force', true);
plot_spec(spec, 'ch', 1, 'logfreq', false, 'saveFig', false,...
    'axh', [])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% psd in states according to timebins
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% analyze
for isession = 1 : nsessions
    basepath = basepaths{isession};
    cd(basepath)
    
    % arrange timebins
    session = v(isession).session;
    timebins = session.general.timebins;
    timepnt = session.general.timepnt;
%     timebins = [1, timepnt; 
%         timepnt, timepnt + 2 * 3600;
%         timepnt + 2 * 3600, timepnt + 4 * 3600;
%         timepnt + 4 * 3600, timepnt + 6 * 3600];

    psdBins = psd_states_timebins('basepath', pwd,...
        'chEeg', [], 'forceA', true, 'graphics', true,...
        'timebins', timebins, 'saveVar', false, 'sstates', [1, 4]);
end

% 2prism
freq = v(1).psdBins.info.freq;
ibin = 2;
istate = 1;
psdState = nan(length(freq), nsessions);
for isession = 1 : nsessions
    psdtmp = squeeze(v(isession).psdBins.psdLfp(ibin, istate, :));
    psdState(:, isession) = psdtmp / sum(psdtmp);
    psdState(:, isession) = psdtmp;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% psd regardless of states
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% analyze
faxis = [0.5 : 0.5 : 120];
clear psd
for isession = 1 : nsessions
    basepath = basepaths{isession};
    [~, basename] = fileparts(basepath);
    cd(basepath)
    
    % arrange winCalc
    session = v(isession).session;
    timebins = session.general.timebins;   
    timepnt = session.general.timepnt;

    timebins = [timepnt, timepnt + 20 * 60];
    
    % load data and calc psd
    ch = v(isession).ss.info.sSig.eegCh;
    for ibin = 1 : size(timebins, 1)
        sig = double(bz_LoadBinary([basename, '.lfp'],...
            'duration', diff(timebins(ibin, :)) + 1,...
            'frequency', 1250, 'nchannels', session.extracellular.nChannels,...
            'start', timebins(ibin, 1), 'channels', ch, 'downsample', 1));
        sig = mean(sig, 2);

        psd(isession, ibin, :) = psd_timebins('sig', sig, 'winCalc', [],...
            'fs', 1250, 'graphics', false, 'faxis', faxis);          
    end
end

ibin = 1;
prism_data = squeeze(psd(:, ibin, :))';
for isession = 1 : nsessions
    prism_data(:, isession) = prism_data(:, isession) / sum(prism_data(:, isession), 1);
end

fh = figure;
plot(faxis, squeeze(mean(psd(:, ibin, :), 1)))
set(gca, 'XScale', 'log')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% time in states
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nbins = 4;
sstates = [1, 4];
totDur = nan(nbins, length(sstates), nsessions);
epLen = cell(length(sstates), nsessions);
for isession = 1 : nsessions
    basepath = basepaths{isession};
    cd(basepath)
    
    [tempDur, epLen(:, isession), timebins] =...
        as_plotZT('nwin', nbins, 'sstates', sstates, 'ss', ss,...
        'graphics', false);
    
    binLen = diff(timebins');
    totDur(:, :, isession) = (tempDur ./ binLen') * 100;
end   

fh = figure;
th = tiledlayout(2, 2, 'TileSpacing', 'Compact');
axh = nexttile;
dataVec = mean(totDur, 3, 'omitnan');
plot([1 : nbins], dataVec(:, istate),...
    'Color', cfg.colors{sstates(istate)}, 'LineWidth', 2)
ylabel('State duration [%]')
ylim([0 100])
ax = gca;
set(ax.YAxis, 'color', cfg.colors{sstates(istate)})


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% to prism
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
