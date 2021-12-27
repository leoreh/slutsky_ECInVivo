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
basepaths{1} = 'G:\RA\hDLX_Gq_WT2\200820_bslDay1';
basepaths{2} = 'G:\RA\hDLX_Gq_Tg\210820_bslDay2Raw2';
basepaths{3} = 'I:\lh86\lh86_210304_070700';                % after op removal
basepaths{4} = 'I:\lh81\lh81_210206_044300';                % saline injection
basepaths{5} = 'K:\Data\lh99\lh99_211219_085802';           % ketamine local
basepaths{6} = 'I:\lh96\lh96_211126_072000';                % ketamine local (s. or)
basepaths{7} = 'F:\Data\Processed\lh96\lh96_211202_070500'; % ketamine local (s.pyr)
basepaths{8} = 'G:\Data\lh93\lh93_210811_102035';           % ketamine local 
basepaths{9} = 'G:\Data\lh95\lh95_210824_083300';           % ketamine local 

% -------------------------------------------------------------------------
% mea
masterpath = 'K:\Data\MEA';
d = dir(masterpath);
d = d([d(:).isdir]);
d = d(~ismember({d(:).name},{'.','..'}));
basenames = {d.name};
nsessions = length(basenames);
for isession = 1 : nsessions
    basepaths{isessions} = fullfile(masterpath, basenames{isession});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% analyze
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

forceA = false;
saveVar = true;

% -------------------------------------------------------------------------
% su
for isession = 1 : nsessions
    
    % load data for analysis
    cd(basepaths{isession})
    [~, basename] = fileparts(basepath);
    spkFile = [basename '.spikes.cellinfo.mat'];
    sessionFile = [basename, '.session.mat'];
    load(spkFile); load(sessionFile)
    
    % CE
    excludeMetrics = [{'monoSynaptic_connections'}, {'deepSuperficial'},...
        {'theta_metrics'}, {'spatial_metrics'}, {'event_metrics'},...
        {'manipulation_metrics'}, {'state_metrics'}, {'psth_metrics'}];
    cell_metrics = ProcessCellMetrics('session', session,...
    'manualAdjustMonoSyn', false, 'summaryFigures', false,...
    'debugMode', false, 'transferFilesFromClusterpath', false,...
    'submitToDatabase', false, 'getWaveformsFromDat', true,...
    'forceReloadSpikes', false, 'showFigures', false,...
    'excludeMetrics', excludeMetrics);
    
    % spike timing metrics
    st = spktimesMetrics('basepath', basepath,...
        'winCalc', [], 'saveVar', saveVar, 'forceA', forceA);
    
    % spike waveform metrics
    swv = spkwvMetrics('basepath', basepath,...
        'saveVar', saveVar, 'forceA', forceA);
    
    % firing rate
    binsize = 60;
    fr = firingRate(spikes.times, 'basepath', basepath,...
        'graphics', false, 'binsize', binsize, 'saveVar', saveVar,...
        'smet', 'GK', 'winBL', winBL, 'winCalc', [0, Inf], 'forceA', forceA);
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
varsFile = ["fr";...
    "mea";...
    "st_metrics";...
    "swv_metrics";...
    "cell_metrics"];
varArray = getSessionVars('dirnames', basenames, 'mousepath', masterpath,...
    'sortDir', false, 'vars', varsFile);

% name of vars for assignment in workspace
vars = ["fr"; "mea"; "st"; "swv"; "cm"];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% concat data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize
tp = []; spkw = []; royer = []; lidor = []; mfr =[]; tau_rise = [];
mizuseki = []; lvr = []; asym = []; hpk = []; rs = []; fs = [];

for isession = 1 : nsessions
    assignVars(varArray, isession, vars)
    
    rs = [rs, selectUnits([], cm, fr, 0, [], [], 'pyr')'];
    fs = [fs, selectUnits([], cm, fr, 0, [], [], 'int')'];
    mfr = [mfr, fr.mfr'];

    asym = [asym, swv.asym];
    hpk = [hpk, swv.hpk];
    tp = [tp, swv.tp];
    spkw = [spkw, swv.spkw];
    
    lvr = [lvr, st.lvr];
    royer = [royer, st.royer];
    lidor = [lidor, st.lidor];
    mizuseki = [mizuseki, st.mizuseki];
    tau_rise = [tau_rise, st.tau_rise];
    
end

mfr = normalize(mfr, 'range', [0.1 1]);
clear units
units(1, :) = logical(rs);
units(2, :) = logical(fs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ---------------------------------------------------------------------
% classification

fh = figure;
subplot(2, 2, 1)
sh = scatter(tp(units(1, :)), royer(units(1, :)),...
    mfr(units(1, :)) * 3000, 'b', '.');
hold on
sh = scatter(tp(units(2, :)), royer(units(2, :)),...
    mfr(units(2, :)) * 3000, 'r', '.');
set(gca, 'yscale', 'log')
xlabel('Trough to Peak [ms]')
ylabel('Burstiness (royer)')
legend({sprintf('RS = %d su', sum(units(1, :))),...
    sprintf('FS = %d su', sum(units(2, :)))})


subplot(2, 2, 2)
sh = scatter(spkw(units(1, :)), tau_rise(units(1, :)),...
    mfr(units(1, :)) * 3000, 'b', '.');
hold on
sh = scatter(spkw(units(2, :)), tau_rise(units(2, :)),...
    mfr(units(2, :)) * 3000, 'r', '.');
set(gca, 'yscale', 'log')
xlabel('Spike Width [ms]')
ylabel('Burstiness (Tau Rise)')

subplot(2, 2, 3)
sh = scatter(asym(units(1, :)), royer(units(1, :)),...
    mfr(units(1, :)) * 3000, 'b', '.');
hold on
sh = scatter(asym(units(2, :)), royer(units(2, :)),...
    mfr(units(2, :)) * 3000, 'r', '.');
set(gca, 'yscale', 'log')
xlabel('Asymmetry [ms]')
ylabel('Burstiness (mizuseki)')

subplot(2, 2, 4)
sh = scatter(hpk(units(1, :)), lvr(units(1, :)),...
    mfr(units(1, :)) * 3000, 'b', '.');
hold on
sh = scatter(hpk(units(2, :)), lvr(units(2, :)),...
    mfr(units(2, :)) * 3000, 'r', '.');
set(gca, 'yscale', 'log')
xlabel('half peak')
ylabel('Irregularity (LvR)')

% save
figname = fullfile(masterpath, 'cellClass');
export_fig(figname, '-jpg', '-transparent', '-r300')

