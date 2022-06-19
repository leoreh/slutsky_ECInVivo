
[basepaths, v, nsessions] = ketInVivo_sessions({'acsf', 'ket', 'ket10'});
[basepaths, v, nsessions] = ketInVivo_sessions({'ket10'});

% params
cfg = as_loadConfig();
nstates = cfg.nstates;
freq = [0.5 : 0.5 : 120];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% per session spec vs. time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

isession = 3;
[expData, xData] = sessions_catVarTime('mname', '',...
    'dataPreset', {'hypnogram', 'spec', 'fr', 'emg_rms'}, 'graphics', true,...
    'basepaths', basepaths(isession), 'xTicksBinsize', 3,...
    'markRecTrans', true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% recalc spec
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% spec
for isession = 1 : nsessions
    
    % file
    basepath = basepaths{isession};
    cd(basepath)
    [~, basename] = fileparts(basepath);
    
    % get eeg channels from sleep signals
    ssfilename = fullfile(basepath, [basename, '.sleep_sig.mat']);
    load(ssfilename, 'info')

    % calc spec
    spec = calc_spec('sig', [], 'fs', 1250, 'graphics', false, 'saveVar', true,...
        'padfft', -1, 'winstep', 5, 'logfreq', true, 'ftarget', [],...
        'ch', [{info.eegCh}], 'force', true);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% psd in states according to timebins
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% analyze
psd = nan(nsessions, 2, length(freq));
for isession = 1 : nsessions
    basepath = basepaths{isession};
    cd(basepath)
    
    % arrange timebins
    session = v(isession).session;
    timebins = session.general.timebins;
    timepnt = session.general.timepnt;
    timebins = [1, timepnt; 
        timepnt, timepnt + 30 * 60];

    psdBins = psd_states_timebins('basepath', pwd,...
        'chEeg', [], 'forceA', true, 'graphics', true,...
        'timebins', timebins, 'saveVar', true, 'sstates', [1]);

    psd(isession, :, :) = squeeze(psdBins.psdLfp);
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
% psd in states according to timebins (legacy)
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

% calc psd
freq = [0.5 : 0.5 : 120];
clear psd
for isession = 1 : nsessions
    basepath = basepaths{isession};
    [~, basename] = fileparts(basepath);
    cd(basepath)

    % arrange timebins
    session = v(isession).session;
    timepnt = session.general.timepnt;
    timebins = [timepnt + 20 * 60, timepnt + 40 * 60];

    % load data
    ch = v(isession).ss.info.sSig.eegCh;
%     ch = [11];
    sig = double(bz_LoadBinary([basename, '.lfp'],...
        'duration', diff(timebins) + 1,...
        'frequency', 1250, 'nchannels', session.extracellular.nChannels,...
        'start', timebins(1), 'channels', ch, 'downsample', 1));
    sig = mean(sig, 2);
    
    % calc psd
    psd(isession, :) = psd_timebins('sig', sig, 'winCalc', [],...
        'fs', 1250, 'graphics', false, 'faxis', freq);
end

% plot psd
fh = figure;
th = tiledlayout(1, 1);
axh = nexttile;
plot(freq, psd)
set(axh, 'XScale', 'log', 'YScale', 'log')
xlabel('Frequency [Hz]')
ylabel('PSD')
legend

% calc bands from psd
bandFreqs = [40, 60];
bandIdx = InIntervals(freq, bandFreqs);
bands = sum(psd(:, bandIdx), 2, 'omitnan');

% to prism
cnt = 1;
prism_data = cell(1, length(grpsz));
for igrp = 1 : length(grpsz)
    prism_data{igrp} = bands(cnt : cnt + grpsz(igrp) - 1);
    cnt = cnt + grpsz(igrp);
end
prism_data = cell2nanmat(prism_data, 2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lfp power in specific bands from spec
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

freq = v(1).spec.freq;
bandFreqs = [0.5, 100; 0.5, 1; 1, 4; 4, 10; 10, 14; 15, 30; 30, 60; 40, 60];
iband = 8;
bands = nan(nsessions, 1);
for isession = 1 : nsessions

    % arrange timebins
    session = v(isession).session;
    timepnt = session.general.timepnt;
    timebins = [timepnt + 1 * 60, timepnt + 20 * 60];
    fs_spec = v(isession).spec.info.winstep;
    timebins = round(timebins / fs_spec);
    tidx = timebins(1) : timebins(2);

    % spec to dB
    powdb = 10 * log10(v(isession).spec.s);
    
    % find frequency indices
    bandIdx = InIntervals(freq, bandFreqs(iband, :));

    bandpow = sum(powdb(:, bandIdx, :), 2);
    bands(isession) = mean(bandpow(tidx), 'omitnan');
end

% to prism
grpsz = [5];
cnt = 1;
prism_data = cell(1, length(grpsz));
for igrp = 1 : length(grpsz)
    prism_data{igrp} = bands(cnt : cnt + grpsz(igrp) - 1, :);
    cnt = cnt + grpsz(igrp);
end
prism_data = cell2nanmat(prism_data, 2);


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
