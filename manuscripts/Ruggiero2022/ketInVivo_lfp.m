
[basepaths, v, nsessions] = ketInVivo_sessions({'acsf', 'ket', 'ket10'});

% params
cfg = as_loadConfig();
nstates = cfg.nstates;
freq = [0.5 : 0.5 : 120];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% per session spec vs. time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

isession = 14;
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
% psd per state for different timebins
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sstates = [1, 4];

% calculate
clear psdBins
for isession = 1 : nsessions
    basepath = basepaths{isession};
    cd(basepath)
    
    % arrange timebins
    session = v(isession).session;
    timepnt = session.general.timepnt;
    timebins = [1, timepnt; timepnt, timepnt + 30 * 60;...
        timepnt, timepnt + 60 * 60; 
        timepnt, timepnt + 120 * 60;
        timepnt, timepnt + 180 * 60];
    timebins = v(isession).session.general.timebins;
%     [timebins, timepnt] = metaInfo_timebins('reqPnt', [], 'nbins', 2);
    timebins = [1, timebins(end)];

    psdBins(isession) = psd_states_timebins('basepath', pwd,...
        'chEeg', [], 'chLfp', [9], 'forceA', true, 'graphics', false,...
        'timebins', timebins, 'saveVar', true, 'sstates', sstates);

end

% cat sessions
psd_cat = catfields([v.psdBins], 'catdef', 'addim', 'force', false);

% organize psd 4D mat [sessions, timebins, freqbands, states]
psd = permute(psd_cat.psdLfp, [4, 1, 3, 2]);    

% -------------------------------------------------------------------------
% organize psd in specific bands and normalize after to before
bandFreqs = [0.5, 1; 1, 4; 4, 10; 30, 60; 60, 100];
nbands = size(bandFreqs, 1);
bands = nan(nsessions, nbands, length(sstates));
for iband = 1 : nbands
    bandIdx = InIntervals(v(1).psdBins(1).info.freq, bandFreqs(iband, :));
    for istate = 1 : length(sstates)
        psd_tmp = sum(psd(:, :, bandIdx, istate), 3);
        bands(:, iband, istate) = psd_tmp(:, 2) ./ psd_tmp(:, 1);
    end
end

% to prism
iband = 2;
istate = 1;
prism_vec = bands(:, iband, istate);
prism_data = vec2nanmat(prism_vec, [6, 5, 5]) * 100;

istate = 2;
ibin = [1, 2];
isession = 5;
x = squeeze(psd(isession, ibin, :, istate))';
movmean(x, 5)
close all
fh = figure;
plot(movmean(x, 5))
set(gca, 'yscale', 'log', 'xscale', 'log')

% stats
setMatlabGraphics(true)
[p, ~, stats] = kruskalwallis(prism_data, {'acsf', 'ket', 'ket10'}, 'on');
c = multcompare(stats, 'Display', 'off');
tbl = array2table(c,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])

% -------------------------------------------------------------------------
% get the raw / norm psd per state
isessions = 12 : 16;
istate = 3;
ibin = 5;
prism_data = squeeze(psd(isessions, ibin, :, istate))';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lfp power in specific bands from spec, with state dependency 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bandFreqs = [0.5, 100; 0.5, 1; 1, 4; 4, 10; 10, 14; 15, 30; 30, 60; 60, 120];
sstates = [1, 4, 5];

nbands = size(bandFreqs, 1);
freq = v(1).spec.freq;
bands = nan(nsessions, 2, nbands, length(sstates));
for isession = 1 : nsessions
    
    % spec to dB
    powdb = 10 * log10(v(isession).spec.s);

    % downsample labels
    fs_spec = v(isession).spec.info.winstep;
    labels = movmedian(v(1).ss.labels, [2 2]);
    labels = labels(1 : fs_spec : end);

    % arrange timebins
    session = v(isession).session;
    timepnt = session.general.timepnt;
    timebins = [1, timepnt; 
        timepnt, timepnt + 30 * 60];
    timebins = round(timebins / fs_spec);
    
    for ibin = 1 : size(timebins, 1)
        for iband = 1 : nbands

            % find frequency indices
            bandIdx = InIntervals(freq, bandFreqs(iband, :));

            % get power in bands
            bandpow = sum(powdb(:, bandIdx, :), 2);

            % find state and time indices, and average psd
            for istate = 1 : length(sstates)
                tidx = intersect(find(labels == istate), timebins(ibin, 1) : timebins(ibin, 2));
                bands(isession, ibin, iband, istate) = mean(bandpow(tidx), 'omitnan');

            end
        end
    end
end

% normalize after to before
bands_norm = nan(nsessions, nbands, length(sstates));
for iband = 1 : nbands
    for istate = 1 : length(sstates)
        bands_norm(:, iband, istate) = bands(:, 2, iband, istate) ./ bands(:, 1, iband, istate);
    end
end

% to prism
iband = 8;
istate = 1;
prism_vec = bands_norm(:, iband, istate);
prism_data = vec2nanmat(prism_vec, [6, 5, 5]);

% stats
setMatlabGraphics(true)
[p, ~, stats] = kruskalwallis(prism_data, {'acsf', 'ket', 'ket10'}, 'on');
c = multcompare(stats, 'Display', 'off');
tbl = array2table(c,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])


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
% legacy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% psd regardless of states
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calc psd
freq = [0.5 : 0.5 : 120];
clear psd
for isession = 12 : 16
    basepath = basepaths{isession};
    [~, basename] = fileparts(basepath);
    cd(basepath)

    % arrange timebins
    session = v(isession).session;
    timepnt = session.general.timepnt;
    timebins = [timepnt + 20 * 60, timepnt + 40 * 60];

    % load data
    ch = v(isession).ss.info.sSig.eegCh;
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