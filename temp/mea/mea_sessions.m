% mea_sessions

% loads all relavent files from multiple sessions (experiments) and does
% stuff

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get all folders in masterpath
masterpath = 'G:\Data\MEA';
d = dir(masterpath);
d = d([d(:).isdir]);
d = d(~ismember({d(:).name},{'.','..'}));
basenames = {d.name};
nsessions = length(basenames);
for isession = 1 : nsessions
    basepaths{isession} = fullfile(masterpath, basenames{isession});
end

% analyze all sessions
for isession = 1 : nsessions
    mea_analyze('basepath', basepaths{isession},...
        'winBL', [0, 120 * 60], 'graphics', false, 'forceA', true)

    % plot fr vs. time
    plot_FRtime_session('basepath', basepaths{isession}, 'grp', [],...
        'frBoundries', [0.01 Inf; 0.01 Inf], 'muFlag', false, 'saveFig', false,...
        'dataType', 'norm')
    
end

% load vars from each session
varsFile = ["fr"; "mea"; "st_metrics"; "swv_metrics"; "cell_metrics"];
varsName = ["fr"; "mea"; "st"; "swv"; "cm"];
varArray = getSessionVars('basepaths', basepaths, 'varsFile', varsFile,...
    'varsName', varsName);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% concat data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize
tp = []; spkw = []; royer = []; lidor = []; mfr =[]; tau_rise = [];
mizuseki = []; lvr = []; asym = []; hpk = []; rs = []; fs = [];

for isession = 1 : nsessions
    
    rs = [rs, selectUnits([], varArray(isession).cm,...
        varArray(isession).fr, 0, [], [], 'pyr')'];
    fs = [fs, selectUnits([], varArray(isession).cm,...
        varArray(isession).fr, 0, [], [], 'int')'];
    mfr = [mfr, varArray(isession).fr.mfr'];

    asym = [asym, varArray(isession).swv.asym];
    hpk = [hpk, varArray(isession).swv.hpk];
    tp = [tp, varArray(isession).swv.tp];
    spkw = [spkw, varArray(isession).swv.spkw];
    
    lvr = [lvr, varArray(isession).st.lvr];
    royer = [royer, varArray(isession).st.royer];
    lidor = [lidor, varArray(isession).st.lidor];
    mizuseki = [mizuseki, varArray(isession).st.mizuseki];
    tau_rise = [tau_rise, varArray(isession).st.tau_rise];
    
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

