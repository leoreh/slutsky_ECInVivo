% wrapper for batch processing

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mname = 'lh133';

varsFile = ["fr"; "sr"; "spikes"; "st_metrics"; "swv_metrics";...
    "cell_metrics"; "sleep_states"; "ripp.mat"; "datInfo"; "session";...
    "units"];
varsName = ["fr"; "sr"; "spikes"; "st"; "swv"; "cm"; "ss"; "ripp";...
    "datInfo"; "session"; "units"];
xlsname = 'D:\Google Drive\PhD\Slutsky\Data Summaries\sessionList.xlsx';
[v, basepaths] = getSessionVars('mname', mname, 'varsFile', varsFile,...
    'varsName', varsName, 'pcond', ["tempflag"], 'ncond', [""],...
    'xlsname', xlsname);
nfiles = length(basepaths);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% analyze data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calData = ss.info.calibrationData;
% for ifile = 1 : nfiles
for ifile = [2, 6]

    % file
    basepath = basepaths{ifile};
    cd(basepath)
    [~, basename] = fileparts(basepath);

    % print progress
    fprintf('sessions_wrapper: working on session %d of %d, %s\n',...
        ifile, nfiles, basename)

    % params
    session = CE_sessionTemplate(pwd, 'viaGUI', false,...
        'forceDef', true, 'forceL', true, 'saveVar', true);
    basepath = session.general.basePath;
    nchans = session.extracellular.nChannels;
    fs = session.extracellular.sr;
    spkgrp = session.extracellular.spikeGroups.channels;
    [~, basename] = fileparts(basepath);


    % calc psd
    ch = [1 : 4];
    sstates = [1, 4];
    ftarget = [0.1 : 0.5 : 100];
    prct = 70;
    flgEmg = true;
    wins = [0, Inf];
    wins = [session.general.timebins(1, 2), Inf]
    v(ifile).psd = psd_states('basepath', basepaths{ifile}, 'sstates', sstates,...
        'ch', ch, 'fs', [], 'wins', wins, 'saveVar', true,...
        'graphics', true, 'forceA', true, 'ftarget', ftarget,...
        'prct', prct, 'flgEmg', flgEmg);


end


% cell_metrics = CellExplorer('basepaths', basepaths);

% concatenate var from different sessions
[expData, xData] = sessions_catVarTime('mname', mname,...
    'dataPreset', {'fr'}, 'graphics', true, 'dataAlt', 1,...
    'basepaths', {}, 'xTicksBinsize', 6, 'markRecTrans', true);

% cell_metrics = CellExplorer('basepath', pwd);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spktimes of single units (rs and fs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%     clear srsu
%     for unitType = [1, 2]
%         spktimes = cell(1, 4);
%         for igrp = 1 : length(spkgrp)
%
%             % get rs and fs indices
%             grpidx = v(ifile).spikes.shankID == igrp;
%             unitidx = v(ifile).units.clean(unitType, :) & grpidx;
%             spktimes{igrp} = sort(vertcat(v(ifile).spikes.times{unitidx}));
%         end
%         srsu(unitType) = calc_fr(spktimes, 'basepath', basepath,...
%             'graphics', false, 'binsize', 60, 'saveVar', false,...
%             'smet', 'none', 'winBL', [0, Inf], 'winCalc', [0, Inf]);
%     end
%     save(fullfile(basepath, [basename, '.srsu.mat']), 'srsu')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% psd across sessions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bands = sessions_psd('lh122', 'flgNormBand', true, 'flgAnalyze', false);









% ---------------- psd from spec

% load vars from each session
varsFile = ["session"; "spec"; "sleep_states"];
varsName = ["session"; "spec"; "ss"];
[v, ~] = getSessionVars('basepaths', basepaths, 'varsFile', varsFile,...
    'varsName', varsName);
nfiles = length(basepaths);

sstates = [1, 4];
clear psd spec
for ifile = 1 : nfiles

    cd(basepaths{ifile})
    basepath = basepaths{ifile};
    [~, basename] = fileparts(basepath);
    nchans = v(ifile).session.extracellular.nChannels;
    ch = 7;
    dur = 6 * 60 * 60;

    sig = mean(double(bz_LoadBinary([basename, '.lfp'], 'duration', dur,...
        'frequency', 1250, 'nchannels', nchans, 'start', 0,...
        'channels', 7, 'downsample', 1)), 2);

    spec(ifile) = calc_spec('basepath', basepaths{ifile},...
        'sig', sig, 'fs', 1250, 'graphics', false,...
        'saveVar', false, 'padfft', -1, 'winstep', 1,...
        'ftarget', [0.1 : 0.5 : 100], 'ch', {[1]},...
        'force', true, 'logfreq', true);

    for istate = 1 : length(sstates)
        stateIdx = v(ifile).ss.labels == sstates(istate);
        stateIdx = stateIdx(1 : dur);
        psd(ifile, istate, :) = mean(spec(ifile).s(stateIdx, :), 1);
    end
end

fh = figure;
th = tiledlayout(1, length(sstates), 'TileSpacing', 'Compact');
alphaIdx = linspace(0.5, 1, nfiles);

for istate = 1 : nstates
    axh = nexttile;

    ydata = squeeze(psd(:, istate, :))';
    ydata = ydata ./ sum(ydata(lim_fAxis, :), 1);

    ph = plot(spec(ifile).freq, ydata, 'LineWidth', 2);
    %     for ifile = 1 : nfiles
    %         ph(ifile).Color(istate) = cfg.colors{sstates(istate)}(istate) - ifile * 0.01;
    %         ph(ifile).Color(4) = alphaIdx(ifile);
    %     end
    set(gca, 'YScale', 'log', 'XScale', 'log')
    title(cfg.names{sstates(istate)})
    xlabel('Frequency [Hz]')
    legend(split(num2str(1 : nfiles)), 'Location', 'Southwest', 'FontSize', 9)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MFR states
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


mname = {'lh123', 'lh122', 'lh126', 'lh129'};

iunit = 1;
istate = 2;

clear fr_states fr_gain states_temp gain_temp
for imouse = 1 : length(mname)
    
    % load data
    varsFile = ["fr"; "units"];
    varsName = ["fr"; "units"];
    xlsname = 'D:\Google Drive\PhD\Slutsky\Data Summaries\sessionList.xlsx';
    [v, basepaths] = getSessionVars('mname', mname{imouse}, 'varsFile', varsFile,...
        'varsName', varsName, 'pcond', ["tempflag"], 'ncond', [""],...
        'xlsname', xlsname);
    nfiles = length(basepaths);

    fr = catfields([v(:).fr], 'catdef', 'cell');

    for ifile = 1 : nfiles
            states_temp{ifile} = fr.states.mfr{ifile}(v(ifile).units.clean(iunit, :), istate)
        gain_temp{ifile} = fr.states.gain{ifile}(4, v(ifile).units.clean(iunit, :))
    end
    fr_states{imouse} = cell2nanmat(states_temp, 2)';
    fr_gain{imouse} = cell2nanmat(gain_temp, 2)';

end

cell2nanmat(fr_states)'
cell2nanmat(fr_gain)'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% burstiness
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load vars from each session
varsFile = ["datInfo"; "session"; "st_metrics"; "units"];
varsName = ["datInfo"; "session"; "st"; "units"];
[v, ~] = getSessionVars('basepaths', basepaths, 'varsFile', varsFile,...
    'varsName', varsName);
nfiles = length(basepaths);

% brst vars to organize and plot
% from avg sub-struct
brstVar = {'freq', 'dur', 'spkprct', 'ibi', 'nspks'};

% fron nbrsts sub-struct

unitType = 'rs';

% initialize
brst = cell2struct(cell(1, length(brstVar)), brstVar, 2);

% concate to cell
for ifile = 1 : nfiles
    su = v(ifile).units.(unitType);

    for ivar = 1 : length(brstVar)
        brst.(brstVar{ivar}){ifile} = v(ifile).st.brst.avg.(brstVar{ivar})(:, su)';
    end

    brst.nbrsts{ifile} = v(ifile).st.brst.nbrsts.freqNorm(:, su)';
end

% organize in mat
for ivar = 1 : length(brstVar)
    brst.(brstVar{ivar}) = cell2nanmat(brst.(brstVar{ivar}), 1);
end
brst.nbrsts = cell2nanmat(brst.nbrsts, 1);

% xlabels
xdata = [-5 * 6 : 6 : 136];

% graphics
setMatlabGraphics(true)
fh = figure;
th = tiledlayout(length(brstVar) + 1, 1, 'TileSpacing', 'Compact');

for ivar = 1 : length(brstVar)
    nexttile
    plot_boxMean('dataMat', brst.(brstVar{ivar}), 'clr', 'k', 'allPnts', true)
    xticklabels(xdata);
    xlabel('Time [h]')
    ylabel(brstVar{ivar})
end

nexttile
plot_boxMean('dataMat', brst.nbrsts, 'clr', 'k', 'allPnts', true)
xticklabels(xdata);
xlabel('Time [h]')
ylabel('Norm. Freq')



