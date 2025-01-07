

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate psd per file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% go over each mouse and analyze all experiment days
grps = [mcu_sessions('wt'), mcu_sessions('mcu')];
flg_bsl = false;

% go over baseline for wt vs. mcu. during experiment, only psdEmg should be
% calculated
grps = {'mcu_bsl'; 'wt_bsl'};
flg_bsl = true;

% state params
ftarget = [1 : 0.5 : 100];

% iterate
for igrp = 1 : length(grps)

    queryStr = grps{igrp};
    basepaths = mcu_sessions(queryStr);
    nfiles = length(basepaths);
    mpath = fileparts(basepaths{1});

    for ifile = 1 : nfiles

        % files
        basepath = basepaths{ifile};
        [~, basename] = fileparts(basepath);
        cd(basepath)
        sigfile = fullfile(basepath, [basename, '.sleep_sig.mat']);
        
        % load lfp signal
        sig = load(sigfile, 'eeg');
        sig = sig.eeg;
        load(sigfile, 'fs');

        % calc psd
        if flg_bsl
            psd = psd_states('basepath', basepath, 'sstates', [1, 4, 5],...
                'sig', sig, 'fs', fs, 'saveVar', true,...
                'graphics', false, 'forceA', true, 'ftarget', ftarget,...
                'flgEmg', false);
        end

        % calc psd accordign to emg state separation
        psdEmg = psd_states('basepath', basepath, 'sstates', [1, 2],...
            'sig', sig, 'fs', fs, 'saveVar', true,...
            'graphics', false, 'forceA', true, 'ftarget', ftarget,...
            'flgEmg', true);
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compare baseline psd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% grab data
[psd_data, psd_cfg] = mcu_psdOrg(2);

% select specific state for analysis
istate = 1;
psd_state = cell(1, 2);
for igrp = 1 : length(psd_data)
    psd_state{igrp} = squeeze(psd_data{igrp}(istate, :, :));
end

% analyze psd difference between groups for this state
[stats, cfg] = mcu_psdCBP(psd_state, psd_cfg);

mcu_psdPlot(stats, psd_state, psd_cfg, 'axh', [], 'type', 'freq',...
    'clr', psd_cfg.clr{istate})


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compare psd across high and low EMG states, across time for WT mice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% grab data
[psd_data, psd_cfg] = mcu_psdOrg(1);
[stats, cfg] = mcu_psdCBPT(psd_data, psd_cfg);

% plot

% open figure
setMatlabGraphics(true)
fh = figure;
set(fh, 'WindowState', 'maximized');
tlayout = [2, ndays];
th = tiledlayout(tlayout(1), tlayout(2));
th.TileSpacing = 'tight';
th.Padding = 'none';
title(th, 'PSD Analysis', 'interpreter', 'none', 'FontSize', 20)
set(fh, 'DefaultAxesFontSize', 16);

% extract params
ndays = size(psd_data{1}, 3);
grps = psd_cfg.snames;
clr = psd_cfg.clr;
faxis = psd_cfg.faxis;
txt_yax = psd_cfg.txt_yax;

axh = nexttile(th, 1, [1, ndays]); cla; hold on
mcu_psdPlot(stats, psd_data, psd_cfg, 'axh', axh, 'type', 'heat', 'flg_log', true)

% plot psd across days
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% psd per state
for iday = 1 : ndays
    axh = nexttile(th, iday + ndays, [1, 1]); cla; hold on
    
    psd_cfg.expDsgn = 2;
    psd_cfg.grps = psd_cfg.snames;
    psd_day = {squeeze(psd_data{1}(:, :, iday))', squeeze(psd_data{2}(:, :, iday))'};
    [stats, ~] = mcu_psdCBPT(psd_day, psd_cfg);
    mcu_psdPlot(stats, psd_day, psd_cfg, 'axh', axh, 'type', 'freq',...
        'flg_log', true)

    title(axh, sprintf('day %d', iday))
end


stats.negclusters(1)




















% normalize to broadband
flg_norm = true;

% can be done using emg or as states
flg_emg = true;

% files
grps = ["wt_bsl"; "mcu_bsl"];
grps = mcu_sessions('wt');

% state params
if flg_emg
    vars = ["psdEmg"];
    sstates = [1, 2];
else
    vars = ["psd"];
    sstates = [1, 4, 5];
end
cfg = as_loadConfig('flgEmg', flg_emg);
nstates = length(sstates);
clr = cfg.colors(sstates);
snames = cfg.names(sstates);

% psd params
if flg_norm
    txt_yax = 'Norm. Power';
else
    txt_yax = 'Power [dB]';
end

% load and organize data
clear psd psdCell
for igrp = 1 : length(grps)
    
    % load
    basepaths = [mcu_sessions(grps{igrp})];
    nfiles = length(basepaths);
    v = basepaths2vars('basepaths', basepaths, 'vars', vars);

    % organize
    psd = catfields([v(:).psd], 'addim', true);    
    psdMat = psd.psd;
    bandsMat = psd.bands.mean;
   
    % normalize
    if flg_norm
        broadPow = sum(psdMat, 2);
        psdMat = psdMat ./ broadPow;
        bandsMat = bandsMat ./ bandsMat(:, :, 1);
    end
    psdCell{igrp} = psdMat;
    bandCell{igrp} = bandsMat;
end

faxis = squeeze(psd.info.faxis(1, :, 1));

% Plot PSD Comparison
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

psd_compare(psdCell, faxis, grps, clr, snames, txt_yax)


% compare 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% reorgnize to cell of states [freq x session x mouse]
stateCell = cell(1, nstates);  

% Get dimensions
nMice = length(psdCell);  % Number of mice
nFreq = size(psdCell{1}, 2);  % Frequency bins
nDays = cellfun(@(x) size(x, 3), psdCell);  % Days per mouse

% Preallocate matrices
for istate = 1 : nstates
    stateCell{istate} = nan(nFreq, max(nDays), nMice);
end

% Rearrange the data
for iMouse = 1:nMice
    % Extract current mouse data (state x freq x day)
    data = psdCell{iMouse};
    
    % For each state, reshape and store
    stateCell{1}(:, 1:size(data, 3), iMouse) = squeeze(data(1, :, :));  % State 1
    stateCell{2}(:, 1:size(data, 3), iMouse) = squeeze(data(2, :, :));  % State 2
end

% run stat analysis
stat = stat_cluPermutation(stateCell, faxis);




figure;
imagesc(stat.time, stat.freq, squeeze(stat.stat));  % Heatmap of t-values
axis xy;
xlabel('Days');
ylabel('Frequency (Hz)');
title('State Differences Across Frequency and Time');
colorbar;

hold on;
contour(stat.time, stat.freq, squeeze(stat.mask), 1, 'k', 'LineWidth', 1.5);  % Significant clusters










% Assume dataCell{1} is the state you want to analyze
data = permute(stateCell{1}, [3, 1, 2]);  % [mouse x freq x day]
nMice = size(data, 1);
nFreq = size(data, 2);
nDays = size(data, 3);




cfg = [];
cfg.method = 'montecarlo';
cfg.statistic = 'depsamplesFunivariate';  % Dependent t-test for repeated measures
cfg.correctm = 'cluster';
cfg.clusteralpha = 0.8;
cfg.tail = 1;  % Two-tailed
cfg.clustertail = 1;
cfg.numrandomization = 1000;

design = zeros(2, nMice);
design(1,:) = 1:nMice;  % Mouse IDs
design(2,:) = ones(1, nMice);  % Condition (all set to 1 initially)


design = zeros(2, nMice * nDays);
design(1,:) = repmat(1:nMice, 1, nDays);  % Mouse IDs (repeat for each day)
design(2,:) = repelem(1:nDays, nMice);    % Day (condition)

design = zeros(2, nMice * nDays);
design(1,:) = repmat(1:nMice, 1, nDays);  % Mouse IDs repeated for each day
design(2,:) = repelem(1:nDays, nMice);    % Condition (day)


cfg.design = design;
cfg.uvar = 1;  % Unit of observation (mouse ID)
cfg.ivar = 2;  % Independent variable (day)


freq_strct_split = cell(1, nDays);
for day = 1:nDays
    freq_strct_split{day}.powspctrm = freq_strct.powspctrm(:,:,:,day);  % Extract data for each day
    freq_strct_split{day}.dimord = 'subj_chan_freq';
    freq_strct_split{day}.freq = faxis;
    freq_strct_split{day}.label = {'dummy_channel'};

end


% Expand to [mouse x 1 x freq x day] (add singleton channel dimension)
freq_strct.powspctrm = reshape(data, nMice, 1, nFreq, nDays);
freq_strct.powspctrm = freq_strct.powspctrm - mean(freq_strct.powspctrm(:,:,:,1), 4);
freq_strct.dimord = 'subj_chan_freq_time';
freq_strct.time = 1:nDays;
freq_strct.freq = faxis;
freq_strct.label = {'dummy_channel'};




stat = ft_freqstatistics(cfg, freq_strct_split{:});


figure;
imagesc(stat.time, stat.freq, squeeze(mean(stat.stat, 1)));  % Average across subjects
xlabel('Days');
ylabel('Frequency (Hz)');
title('Cluster Permutation Test Across Days (depsamplesFunivariate)');
colorbar;












%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compare psd from emg and as
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% normalize to broadband
flg_norm = true;

% files
grps = ["wt_bsl"; "mcu_bsl"];


% state params
if flg_emg
    sstates = [1, 2];
else
    vars = ["psd"];
    sstates = [1, 4, 5];
end
cfg = as_loadConfig('flgEmg', flg_emg);
nstates = length(sstates);
clr = cfg.colors(sstates);
snames = cfg.names(sstates);

% psd params
if flg_norm
    txt_yax = 'Norm. Power';
else
    txt_yax = 'Power [dB]';
end





% load and organize data
clear psd psdCell
    
% load
basepaths = [mcu_sessions('wt_bsl'), mcu_sessions('mcu_bsl')];
nfiles = length(basepaths);

% as
vars = ["psd"];
v = basepaths2vars('basepaths', basepaths, 'vars', vars);
psd = catfields([v(:).psd], 'addim', true);
psdMat = psd.psd([1, 2], :, :);
if flg_norm
    broadPow = sum(psdMat, 2);
    psdMat = psdMat ./ broadPow;
end
psdCell{1} = psdMat;

% emg
vars = ["psdEmg"];
v = basepaths2vars('basepaths', basepaths, 'vars', vars);
psd = catfields([v(:).psd], 'addim', true);
psdMat = psd.psd;
if flg_norm
    broadPow = sum(psdMat, 2);
    psdMat = psdMat ./ broadPow;
end
psdCell{2} = psdMat;

faxis = squeeze(v(1).psd.info.faxis(1, :, 1));

% Plot PSD Comparison
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

psd_compare(psdCell, faxis, grps, clr, snames, txt_yax)








% To do
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1 % investigate clu permute with time
% 1 % compare high to low emg across time 
% 1 % compare mcu to wt across time for each state

% plot heat matrix of t-values with contour of significant






% 
% % Plot Bands Comparison
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% bands_idx = [2 :  6];
% nbands = length(bands_idx);
% 
% % open figure
% setMatlabGraphics(true)
% fh = figure;
% set(fh, 'WindowState', 'maximized');
% tlayout = [nstates, nbands];
% th = tiledlayout(tlayout(1), tlayout(2));
% th.TileSpacing = 'tight';
% th.Padding = 'none';
% title(th, 'Bands Baseline', 'interpreter', 'none', 'FontSize', 20)
% set(fh, 'DefaultAxesFontSize', 16);
% 
% % psd per state
% for istate = 1 : nstates
%     for iband = 1 : nbands
% 
%         % organize data
%         clear dataMat
%         for igrp = 1 : length(grps)
%             dataMat{igrp} = squeeze(bandCell{igrp}(istate, bands_idx(iband), :));
%         end
%         dataMat = cell2padmat(dataMat, 2);
% 
%         % plot
%         axh = nexttile(th, nbands * (istate - 1) + iband, [1, 1]); cla; hold on
%         plot_boxMean('axh', axh, 'dataMat', dataMat, 'allPnt', false,...
%         'plotType', 'bar');
%         title(axh, psd.bands.info.bandNames{1, bands_idx(iband), 1})
%         xticklabels({'WT', 'MCU-KO'})
%     end
% end










% 
% % get psd across sessions for one mouse
% [bands, powdb] = sessions_psd('mname', mname, 'flgNormBand', true, 'flgDb', false,...
%     'flgNormTime', false, 'flgEmg', true, 'idxBsl', [1 : 2], 'saveFig', false);
