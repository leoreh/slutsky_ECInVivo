
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data base
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% aCSF
basepaths = [
    {'F:\Data\lh93\lh93_210811_102035'},...
    {'F:\Data\lh95\lh95_210824_083300'},...
    {'F:\Data\lh96\lh96_211201_070100'},...
    {'F:\Data\lh99\lh99_211218_090630'},...
    {'F:\Data\lh100\lh100_220405_100406'},...
    {'F:\Data\lh107\lh107_220509_095738'},...
    ];
%     {'F:\Data\Processed\lh96\lh96_211124_073800'},...
%     {'F:\Data\lh93\lh93_210811_102035'},...
%     {'F:\Data\lh99\lh99_211218_090630'},...

% local ket
basepaths = [
    {'F:\Data\lh93\lh93_210813_110609'},...
    {'F:\Data\lh95\lh95_210825_080400'},...
    {'F:\Data\lh96\lh96_211204_084200'},...
    {'F:\Data\lh100\lh100_220403_100052'},...
    {'F:\Data\lh107\lh107_220501_102641'},...
    ];
% {'F:\Data\lh99\lh99_220119_090035'},...
%     {'F:\Data\lh96\lh96_211126_072000'},...
%     {'F:\Data\lh96\lh96_211202_070500'},...

% ip ket 10 mg/kg
basepaths = [
    {'F:\Data\lh96\lh96_211206_070400'},...
    {'F:\Data\lh98\lh98_211224_084528'},...
    {'F:\Data\lh106\lh106_220512_102302'},...
    {'F:\Data\lh107\lh107_220512_102302'},...
    ];
% lh99_211224_084528
% F:\Data\lh81\lh81_210204_190000   inj 0zt

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reanalyze something
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load vars from each session
forceL = true;
varsFile = ["fr"; "fr_bins"; "spikes"; "datInfo"; "session"; "units"];
varsName = ["fr"; "frBins"; "spikes"; "datInfo"; "session"; "units"];
vload = getSessionVars('basepaths', basepaths, 'varsFile', varsFile,...
    'varsName', varsName);
nsessions = length(basepaths);

cell_metrics = CellExplorer('basepaths', basepaths);

% itterate
for isession = 1 : nsessions
    basepath = basepaths{isession};
    cd(basepath)
    [mousename, basename] = fileparts(basepath);
    [~, mousename] = fileparts(mousename);

    % number of units per spike group
    plot_nunits_session('basepath', basepath, 'frBoundries', [])

    % select specific units
    units = selectUnits('basepath', basepath, 'grp', [2 : 7], 'saveVar', true,...
        'forceA', true, 'frBoundries', [0 Inf; 0 Inf],...
        'spikes', vload(isession).spikes);
    units = vload(isession).units;
    unitsClean = units.idx & units.gini' & units.stable' & units.mfrBL' &...
        units.su' & units.cnt';
    
  
end

% params
cfg = as_loadConfig();
nstates = cfg.nstates;
unitChar = {'RS', 'FS'};
unitClr = {'b', 'r'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% concate sessions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load vars from each session
varsFile = ["session"; "units"; "fr_bins"; "fr"; "sleep_states"];
varsName = ["session"; "units"; "frBins"; "fr"; "ss"];
v = getSessionVars('basepaths', basepaths, 'varsFile', varsFile,...
    'varsName', varsName);
nsessions = length(basepaths);

% select states
sstates = [1, 4];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics - psd in states according to timebins
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% analyze
for isession = 1 : nsessions
    cd(basepaths{isession});
    session = v(isession).session;
    [timebins, timepnt] = metaInfo_timebins('reqPnt', 5.5 * 60 * 60,...
        'nbins', 8);
    psdBins = psd_states_timebins('basepath', pwd,...
        'chEeg', [], 'forceA', true, 'graphics', true,...
        'timebins', timebins, 'saveVar', true, 'sstates', [1, 4, 5]);
end

% load
varsFile = ["psdBins"; "session"];
varsName = ["psdBins"; "session"];
v = getSessionVars('basepaths', basepaths, 'varsFile', varsFile,...
    'varsName', varsName);

% 2prism
freq = psdBins.info.freq;
ibin = 2;
istate = 3;
psdState = nan(length(freq), nsessions);
for isession = 1 : nsessions
    psdtmp = squeeze(v(isession).psdBins.psdLfp(ibin, istate, :));
    psdState(:, isession) = psdtmp / sum(psdtmp);
    psdState(:, isession) = psdtmp;
end

% time in states
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

% firing rate vs. time across all units, aligned to point of injection
units = catfields([v.units], 'catdef', 'long', 'force', false);
unitsClean = units.idx & units.gini' & units.stable' & units.mfrBL' &...
    units.su' & units.cnt';

% firing rate vs. time across all units, aligned to point of injection. can
% extract strd or normalized data
[frMat, timeIdx] = alignFR2pnt('basepaths', basepaths, 'dataType', 'strd');

% grab fr
unitNo = 1;
prism_data = frMat(unitsClean(unitNo, :), :)';
prism_tstamps = [1 : length(frMat)]' / 60;
prism_injTime = max(timeIdx) / 60;

% list of units per session
prism_rs = []; prism_fs = [];
for isession = 1 : nsessions
    units = v(isession).units;
    unitsClean = units.idx & units.gini' & units.stable' & units.mfrBL' &...
    units.su' & units.cnt';
    prism_rs = [prism_rs; ones(sum(unitsClean(1, :)), 1) * isession];
    prism_fs = [prism_fs; ones(sum(unitsClean(2, :)), 1) * isession];
end

% select units based on mfr
fr = catfields([v.fr], 'catdef', 'long', 'force', false);
units_mfr = fr.mfr < 3;
unitsClean = unitsClean & units_mfr';

% single session
dataType = 'norm';
data = fr.(dataType)(units.idx(unitNo, :), :)';
data(~isfinite(data)) = nan;
excludeIdx = ~units.gini' | ~units.mfrBL';
prismIdx = excludeIdx(units.idx(unitNo, :));
tstamps = fr.tstamps / 60 / 60;