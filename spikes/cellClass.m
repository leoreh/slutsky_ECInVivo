% cellClass

% wrapper. gets multiple basepaths (manually or all paths in a master
% folder), anaylzes each session for spike times, waveforms metrices, and
% firing rate. loads and concatenates the metrices. plots the metrices with the
% expectation that two clusters (rs and fs) will be easily visualized. 
% next step to is to apply a classifying algorithm (e.g. GMM) on selected
% matrices (e.g. tp, spkw, and tau_rise)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% session paths that will serve as the data base
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -------------------------------------------------------------------------
% single units 
basepaths{1} = 'F:\Data\Colleagues\RA\tg4_210730_155700';
basepaths{2} = 'F:\Data\Colleagues\RA\hDLX_Gq_WT5\040221_0655_24hr';
basepaths{3} = 'I:\lh86\lh86_210304_070700';                % after op removal
basepaths{4} = 'I:\lh81\lh81_210206_044300';                % saline injection
basepaths{5} = 'K:\Data\lh99\lh99_211219_085802';           % ketamine local
basepaths{6} = 'I:\lh96\lh96_211126_072000';                % ketamine local (s. or)
basepaths{7} = 'F:\Data\Processed\lh96\lh96_211202_070500'; % ketamine local (s.pyr)
basepaths{8} = 'G:\Data\lh93\12hr\lh93_210811_102035';      % ketamine local 
basepaths{9} = 'K:\Data\lh95\12hr\lh95_210824_083300';      % ketamine local 

% -------------------------------------------------------------------------
% mea
masterpath = 'K:\Data\MEA';
d = dir(masterpath);
d = d([d(:).isdir]);
d = d(~ismember({d(:).name},{'.','..'}));
basenames = {d.name};
nsessions = length(basenames);
for isession = 1 : nsessions
    basepaths{isession} = fullfile(masterpath, basenames{isession});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% analyze
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

forceA = true;
saveVar = true;
nsessions = length(basepaths);

% -------------------------------------------------------------------------
% su
for isession = 1 : nsessions
    
    % load data for analysis
    basepath = basepaths{isession};
    cd(basepath)
    [~, basename] = fileparts(basepath);
    spkFile = [basename '.spikes.cellinfo.mat'];
    sessionFile = [basename, '.session.mat'];
    load(spkFile); load(sessionFile)
    
    % spike timing metrics
    st = spktimesMetrics('basepath', basepath,...
        'winCalc', [], 'saveVar', saveVar, 'forceA', false);
    
    % spike waveform metrics
    swv = spkwvMetrics('basepath', basepath,...
        'saveVar', saveVar, 'forceA', true);
    
    % firing rate
    binsize = 60;
    winBL = [0 120 * 60];
    fr = firingRate(spikes.times, 'basepath', basepath,...
        'graphics', true, 'binsize', binsize, 'saveVar', saveVar,...
        'smet', 'GK', 'winBL', winBL, 'winCalc', [0, Inf], 'forceA', false);
end

% -------------------------------------------------------------------------
% mea
for isession = 1 : nsessions
    mea_analyze('basepath', fullfile(masterpath, basenames{isession}),...
        'winBL', [0, 120 * 60], 'graphics', false)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load vars from each session
varsFile = ["fr"; "mea"; "st_metrics"; "swv_metrics"; "cell_metrics"];
varsName = ["fr"; "mea"; "st"; "swv"; "cm"];
varArray = getSessionVars('basepaths', basepaths, 'varsFile', varsFile,...
    'varsName', varsName);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% concat data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize
royer = []; lidor = []; tau_rise = []; mizuseki = []; lvr = []; 
tp = []; asym = []; hpk = []; ampP = []; ampTp = []; tslope = []; spkw = []; rtau = [];
mfr =[]; pyr = []; int = [];

for isession = 1 : nsessions
    
    pyr = [pyr, selectUnits([], varArray(isession).cm,...
        varArray(isession).fr, 0, [], [], 'pyr')'];
    int = [int, selectUnits([], varArray(isession).cm,...
        varArray(isession).fr, 0, [], [], 'int')'];
    mfr = [mfr, varArray(isession).fr.mfr'];
    
    asym = [asym, varArray(isession).swv.asym];
    hpk = [hpk, varArray(isession).swv.hpk];
    tp = [tp, varArray(isession).swv.tp];
    spkw = [spkw, varArray(isession).swv.spkw];
    rtau = [rtau, varArray(isession).swv.spkw];
    tslope = [tslope, varArray(isession).swv.spkw];
    ampP = [ampP, varArray(isession).swv.spkw];
    ampTp = [ampTp, varArray(isession).swv.spkw];

    lvr = [lvr, varArray(isession).st.lvr];
    royer = [royer, varArray(isession).st.royer];
    lidor = [lidor, varArray(isession).st.lidor];
    mizuseki = [mizuseki, varArray(isession).st.mizuseki];
    tau_rise = [tau_rise, varArray(isession).st.tau_rise];
    
end

mfr = normalize(mfr, 'range', [0.1 1]);
clear units
units(1, :) = logical(pyr);
units(2, :) = logical(int);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cell_metrics = CellExplorer('basepaths', basepaths);

% ---------------------------------------------------------------------
% classification

fh = figure;
subplot(3, 2, 1)
sh = scatter(tp(units(1, :)), royer(units(1, :)),...
    mfr(units(1, :)) * 3000, 'b', '.');
hold on
sh = scatter(tp(units(2, :)), royer(units(2, :)),...
    mfr(units(2, :)) * 3000, 'r', '.');
set(gca, 'yscale', 'log')
xlabel('Trough to Peak [ms]')
ylabel('Burstiness (royer)')
legend({sprintf('RS = %d su', sum(units(1, :))),...
    sprintf('FS = %d su', sum(units(2, :)))}, 'Location', 'southeast')


subplot(3, 2, 2)
sh = scatter(spkw(units(1, :)), tau_rise(units(1, :)),...
    mfr(units(1, :)) * 3000, 'b', '.');
hold on
sh = scatter(spkw(units(2, :)), tau_rise(units(2, :)),...
    mfr(units(2, :)) * 3000, 'r', '.');
set(gca, 'yscale', 'log')
xlabel('Spike Width [ms]')
ylabel('Burstiness (Tau Rise)')

subplot(3, 2, 3)
sh = scatter(asym(units(1, :)), royer(units(1, :)),...
    mfr(units(1, :)) * 3000, 'b', '.');
hold on
sh = scatter(asym(units(2, :)), royer(units(2, :)),...
    mfr(units(2, :)) * 3000, 'r', '.');
set(gca, 'yscale', 'log')
xlabel('Asymmetry [ms]')
ylabel('Burstiness (mizuseki)')

subplot(3, 2, 4)
sh = scatter(hpk(units(1, :)), lvr(units(1, :)),...
    mfr(units(1, :)) * 3000, 'b', '.');
hold on
sh = scatter(hpk(units(2, :)), lvr(units(2, :)),...
    mfr(units(2, :)) * 3000, 'r', '.');
set(gca, 'yscale', 'log')
xlabel('half peak')
ylabel('Irregularity (LvR)')

subplot(3, 2, 5)
sh = scatter(tslope(units(1, :)), tp(units(1, :)),...
    mfr(units(1, :)) * 3000, 'b', '.');
hold on
sh = scatter(tslope(units(2, :)), tp(units(2, :)),...
    mfr(units(2, :)) * 3000, 'r', '.');
set(gca, 'yscale', 'log')
xlabel('tail slope')
ylabel('trough to peak')

subplot(3, 2, 6)
sh = scatter(tp(units(1, :)), rtau(units(1, :)),...
    mfr(units(1, :)) * 3000, 'b', '.');
hold on
sh = scatter(tp(units(2, :)), rtau(units(2, :)),...
    mfr(units(2, :)) * 3000, 'r', '.');
set(gca, 'yscale', 'log')
xlabel('trough to peak')
ylabel('Repolarization time')

% save
figname = fullfile(masterpath, 'cellClass');
export_fig(figname, '-jpg', '-transparent', '-r300')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% examine waveform normalization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% check the trough aligment of the waveforms

% compare waveform normalization 
basepath = 'K:\Data\lh95\12hr\lh95_210824_083300';
[~, basename] = fileparts(basepath);
load(fullfile(basepath, [basename, '.swv_raw.mat']))
load(fullfile(basepath, [basename, '.cell_metrics.cellinfo.mat']))

fs = 24414.0625;

nunits = 10;
for iunit = 1 : nunits
    
    nspks = size(swv_raw{iunit}, 2);
    for ispk = 1 : nspks
        v = swv_raw{iunit}(:, ispk);
        wv_l2(:, ispk) = v / vecnorm(v, 2, 1);
        wv_norm(:, ispk) = v / abs(min(v));
    end
    
    
    fh = figure;
    x_val = [1 : length(v)] / fs * 1000;

    subplot(2, 2, 1)
    wv_mean = mean(swv_raw{iunit}, 2)';
    wv_std = std(swv_raw{iunit}, [], 2)';
    plot(x_val, wv_mean);
    patch([x_val, flip(x_val)], [wv_mean + wv_std, flip(wv_mean - wv_std)],...
        'k', 'EdgeColor', 'none', 'FaceAlpha', .2, 'HitTest', 'off')
    xlabel('Time [ms]')
    ylabel('Voltage')
    
    subplot(2, 2, 2)
    wv_mean = mean(wv_l2, 2)';
    wv_std = std(wv_l2, [], 2)';
    plot(x_val, wv_mean);
    patch([x_val, flip(x_val)], [wv_mean + wv_std, flip(wv_mean - wv_std)],...
        'k', 'EdgeColor', 'none', 'FaceAlpha', .2, 'HitTest', 'off')
    xlabel('Time [ms]')
    ylabel('L2')
    
    subplot(2, 2, 3)
    wv_mean = cell_metrics.waveforms.raw{iunit};
    wv_std = cell_metrics.waveforms.raw_std{iunit};
    plot(x_val, wv_mean);
    patch([x_val, flip(x_val)], [wv_mean + wv_std, flip(wv_mean - wv_std)],...
        'k', 'EdgeColor', 'none', 'FaceAlpha', .2, 'HitTest', 'off')
    xlabel('Time [ms]')
    ylabel('Normalized')
    subtitle('CE raw')
    
    subplot(2, 2, 4)
    wv_mean_raw = mean(swv_raw{iunit}, 2)';
    wv_mean_l2 = mean(wv_l2, 2)';
    wv_l2_avg = wv_mean_raw ./ vecnorm(wv_mean_raw, 2, 2);
    wv_mean_norm = mean(wv_norm, 2)';

    wv_mean_norm = mean(wv_norm, 2)';
    plot(x_val, wv_mean_raw);
    yyaxis right
    hold on
    plot(x_val, wv_mean_l2);
    plot(x_val, wv_l2_avg);
    plot(x_val, wv_mean_norm);
    legend({'Raw', 'L2', 'L2_avg', 'norm'}, 'Location', 'southeast')
    xlabel('Time [ms]')
    ylabel('Voltage')
    axis tight
    
end


