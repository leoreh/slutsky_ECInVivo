

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% analysis per mouse
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mname = 'lh96';

varsFile = ["fr"; "datInfo"; "session"; "units"];
varsName = ["fr"; "datInfo"; "session"; "units"];
xlsname = 'D:\Google Drive\PhD\Slutsky\Data Summaries\sessionList.xlsx';
[v, basepaths] = getSessionVars('mname', mname, 'varsFile', varsFile,...
    'varsName', varsName, 'pcond', ["tempflag"], 'ncond', [""],...
    'xlsname', xlsname);
nfiles = length(basepaths);


% fr in time bins 
for ifile = 1 : nfiles

    % file
    basepath = basepaths{ifile};
    cd(basepath)
    [~, basename] = fileparts(basepath);

    % print progress
    fprintf('sessions_wrapper: working on session %d of %d, %s\n',...
        ifile, nfiles, basename)

    % add timebins to datInfo
    [timebins, timepnt] = metaInfo_timebins('reqPnt', [], 'nbins', 4);
    timebins / 60 / 60
    
    plot_FRtime_session('basepath', basepath)

    % mfr by states in time bins
    frBins = fr_timebins('basepath', pwd, 'forceA', true, 'graphics', true,...
        'timebins', timebins, 'saveVar', true, 'sstates', [1, 4, 5]);

end


% other stuff
for ifile = 1 : nfiles

    % file
    basepath = basepaths{ifile};
    cd(basepath)
    [~, basename] = fileparts(basepath);

    % print progress
    fprintf('sessions_wrapper: working on session %d of %d, %s\n',...
        ifile, nfiles, basename)

    % select specific units
    load([basename, '.spikes.cellinfo.mat'])
    units = selectUnits('basepath', pwd, 'grp', [1 : 4], 'saveVar', true,...
        'forceA', true, 'frBoundries', [0.0 Inf; 0.0 Inf],...
        'spikes', spikes, 'altClean', 2);

    fr = calc_fr(spikes.times, 'basepath', basepath,...
        'graphics', true, 'binsize', 60, 'saveVar', true, 'forceA', true,...
        'smet', 'none', 'winBL', [0 Inf], 'winCalc', [0, Inf]);

end

cell_metrics = CellExplorer('basepaths', basepaths);


% concatenate var from different sessions
[expData, xData] = sessions_catVarTime('mname', mname,...
    'dataPreset', {'sr', 'fr', 'spec'}, 'graphics', true, 'dataAlt', 1,...
    'basepaths', {}, 'xTicksBinsize', 6, 'markRecTrans', true);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get MFR in time bins per unit across all mice, per state or irrespective
% of states
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mname = {'lh122'; 'lh123'; 'lh126'; 'lh129'; 'lh130'};

sstates = [1, 4];
mfrcat = cell(length(sstates) + 1, 2);
nbins = 4;

for imouse = 1 : length(mname)

    % get basepaths
    queryStr = [mname{imouse}, '_mk801'];
    basepaths = mk801_chronic_sessions(queryStr);

    % reload data
    varsFile = ["fr"; "fr_bins"; "datInfo"; "session"; "units"];
    varsName = ["fr"; "frBins"; "datInfo"; "session"; "units"];
    xlsname = 'D:\Google Drive\PhD\Slutsky\Data Summaries\sessionList.xlsx';
    [v, basepaths] = getSessionVars('basepaths', basepaths, 'varsFile', varsFile,...
        'varsName', varsName, 'pcond', ["tempflag"], 'ncond', [""],...
        'xlsname', xlsname);
    nfiles = length(basepaths);

    % organize in cell array
    cnt = 1;
    mfr = cell(nfiles * nbins, length(sstates) + 1, 2);
    for ifile = 1 : nfiles
        for ibin = 1 : nbins
            for iunit = 1 : 2

                unitIdx = v(ifile).units.clean(iunit, :);

                for istate = 1 : length(sstates)                                      
                    mfr{cnt, istate, iunit} = v(ifile).frBins(ibin).states.mfr(unitIdx, sstates(istate));
                end
                mfr{cnt, length(sstates) + 1, iunit} = v(ifile).frBins(ibin).mfr(unitIdx);
            end
            cnt = cnt + 1;
        end
    end

    % reorganize for prism
    for istate = 1 : length(sstates) + 1
        for iunit = 1 : 2
            data = cell2nanmat(squeeze(mfr(:, istate, iunit)), 2);
            mfrcat{istate, iunit} = [mfrcat{istate, iunit}; data];
        end
    end
end
timebins = [-30 : 6 : 84];
tIdx = find(ismember(timebins, [-12, 12, 42, 54]));

% to prism: mfr in bins 
iunit = 1;
istate = 3;
mfrcat{istate, iunit}

% graphics
fh = figure;
th = tiledlayout(2, 2);
xi = linspace(-5, 2, 100);
lgdTxt = {'BSL', 'Acute', 'Chronic', 'WASH'};
for ibin = 1 : 4
    axh = nexttile(th, ibin, [1, 1]); cla; hold on
    data = log10(mfrcat{istate, iunit}(:, tIdx(ibin)));
    histogram(data, 10, 'Normalization', 'pdf')
    [f(ibin, :), xi] = ksdensity(data, xi, 'Bandwidth', 0.4);
    plot(xi, f(ibin, :), 'LineWidth', 2);
    title(axh, lgdTxt{ibin})
    xlabel('log MFR (Hz)')
    ylabel('No. Units (Hz)')
    ylim([0 1])
end


plot_boxMean(mfrcat{istate, iunit})
histogram(mfr)




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% distribution of firing rate bins from each timebin per unit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% also estimates the disribution of firing rate bins for each unit

mname = {'lh123'; 'lh126'; 'lh129'; 'lh130'};

sstates = [1, 4];
mfrcat = cell(length(sstates) + 1, 2);
nbins = 4;
xi = linspace(-3.5, 1.5, 100);
dists = cell(nfiles * nbins, 2);

for imouse = 1 : length(mname)

    % get basepaths
    queryStr = [mname{imouse}, '_mk801'];
    basepaths = mk801_chronic_sessions(queryStr);

    % reload data
    varsFile = ["fr"; "fr_bins"; "datInfo"; "session"; "units"];
    varsName = ["fr"; "frBins"; "datInfo"; "session"; "units"];
    xlsname = 'D:\Google Drive\PhD\Slutsky\Data Summaries\sessionList.xlsx';
    [v, basepaths] = getSessionVars('basepaths', basepaths, 'varsFile', varsFile,...
        'varsName', varsName, 'pcond', ["tempflag"], 'ncond', [""],...
        'xlsname', xlsname);
    nfiles = length(basepaths);

    % organize in cell array
    cnt = 1;
    mfr = cell(nfiles * nbins, length(sstates) + 1, 2);
    for ifile = 1 : nfiles
        for ibin = 1 : nbins
            for iunit = 1 : 2

                unitIdx = v(ifile).units.clean(iunit, :);

                for istate = 1 : length(sstates)                                      
                    data = v(ifile).frBins(ibin).states.fr{istate}(unitIdx, :);
                    mfr{cnt, istate, iunit} = data(:);
                end
                data = v(ifile).frBins(ibin).strd(unitIdx, :);
                mfr{cnt, length(sstates) + 1, iunit} = data(:);

                data = log10(data);
                tmp = nan(size(data, 1), length(xi));
                for icell = 1 : size(data, 1)
                    tmp(icell, :) = ksdensity(data(icell, :), xi, 'Bandwidth', 0.3);
                end
                dists{cnt, iunit} = [dists{cnt, iunit}; tmp];
            end
            cnt = cnt + 1;
        end
    end

    % reorganize for prism
    for istate = 1 : length(sstates) + 1
        for iunit = 1 : 2
            data = cell2nanmat(squeeze(mfr(:, istate, iunit)), 2);
            mfrcat{istate, iunit} = [mfrcat{istate, iunit}; data];
        end
    end
end
timebins = [-30 : 6 : 84];
tIdx = find(ismember(timebins, [-12, 12, 42, 54]));

% to prism: mfr in bins 
iunit = 1;
istate = 1;
mfrcat{istate, iunit};

% graphics
fh = figure;
th = tiledlayout(2, 2);
lgdTxt = {'BSL', 'Acute', 'Chronic', 'WASH'};
for ibin = 1 : 4
    axh = nexttile(th, ibin, [1, 1]); cla; hold on
    data = dists{tIdx(ibin), iunit};
    plot(xi, data, 'LineWidth', 1);
    
    plot(xi, mean(data, 1, 'omitnan'), 'k', 'LineWidth', 3)
    
    title(axh, lgdTxt{ibin})
    xlabel('log MFR (Hz)')
    ylabel('No. Units (Hz)')
%     ylim([0 1])
end

% graphics
fh = figure;
th = tiledlayout(2, 2);
lgdTxt = {'BSL', 'Acute', 'Chronic', 'WASH'};
for ibin = 1 : 4
    axh = nexttile(th, ibin, [1, 1]); cla; hold on
    data = log10(mfrcat{istate, iunit}(:, tIdx(ibin)));
    histogram(data, 'Normalization', 'pdf')
    [f(ibin, :), xi] = ksdensity(data, xi, 'Bandwidth', 0.4);
    plot(xi, f(ibin, :), 'LineWidth', 2);
    title(axh, lgdTxt{ibin})
    xlabel('log MFR (Hz)')
    ylabel('No. Units (Hz)')
    ylim([0 1])
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% re-calc firing rates in smaller time bins
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mname = {'lh123'; 'lh126'; 'lh129'; 'lh130'};

iunit =  1;
binsize = 1;      % [s]
tstamps = [0 : binsize : (5 * 24 * 60 * 60) - 0.1] - (30 * 60 * 60);
recLen = (24 * 60 * 60) / binsize;

catMat = cell(length(mname), 1);
for imouse = 1 : length(mname)

    % get basepaths
    queryStr = [mname{imouse}, '_mk801'];
    basepaths = mk801_chronic_sessions(queryStr);

    % reload data
    varsFile = ["spikes"; "datInfo"; "session"; "units"];
    varsName = ["spikes"; "datInfo"; "session"; "units"];
    xlsname = 'D:\Google Drive\PhD\Slutsky\Data Summaries\sessionList.xlsx';
    [v, basepaths] = getSessionVars('basepaths', basepaths, 'varsFile', varsFile,...
        'varsName', varsName, 'pcond', ["tempflag"], 'ncond', [""],...
        'xlsname', xlsname);
    nfiles = length(basepaths);

    for ifile = 1 : nfiles
        basepath = char(basepaths{ifile});
        
        unitIdx = v(ifile).units.clean(iunit, :);

        % firing rate
        if isfield(v(ifile).session.general, 'timepnt')
            timepnt = v(ifile).session.general.timepnt;
        else
            timepnt = Inf;
        end
        winBL = [0 timepnt];
        fr = calc_fr(v(ifile).spikes.times, 'basepath', basepath,...
            'graphics', false, 'binsize', binsize, 'saveVar', false, 'forceA', true,...
            'smet', 'none', 'winBL', winBL, 'winCalc', [0, Inf]);
        
        data = fr.strd(unitIdx, :);
        [m, n] = size(data);
        if n < recLen
            segLen = recLen - n;
            lastSegIdx = max(1, n - segLen + 1) : n;
            numRepeats = ceil(segLen / numel(lastSegIdx));
            repeatedSegments = repmat(data(:, lastSegIdx), 1, numRepeats);
            % Concatenate the original matrix with the needed portion of the repeated segments
            data = [data, repeatedSegments(:, 1 : segLen)];
        elseif n > recLen
            data = data(:, 1 : recLen);
        end
        frMat{ifile} = data;

    end
    
    % organize in mat
    catMat{imouse} = cell2nanmat(frMat, 2);

end

[newMat, cpad] = cell2nanmat(catMat, 1);


fh = figure;
plot(tstamps / 60 / 60, mean(newMat, 1, 'omitnan'))

prismIdx = round(linspace(1, size(newMat, 2), 1000));
prismMat = newMat(:, prismIdx)';
[tstamps(prismIdx) / 60 / 60]';

mean(newMat(:, prismIdx), 1, 'omitnan')';
std(newMat(:, prismIdx), [], 1, 'omitnan')';
repmat(size(newMat, 1), 1, size(newMat, 2))';
[tstamps / 60 / 60]';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MFR in states across mice per unit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CHECK SESSION_LIST.XLSX

sstates = [1 : 6];
mname = {'lh122'; 'lh123'; 'lh126'; 'lh129'; 'lh130'};
clear mfr
for imouse = 1 : length(mname)
    
    queryStr = [mname{imouse}, '_mk801'];
    basepaths = mk801_chronic_sessions(queryStr);
    
    varsFile = ["fr"; "datInfo"; "session"; "units"];
    varsName = ["fr"; "datInfo"; "session"; "units"];
    xlsname = 'D:\Google Drive\PhD\Slutsky\Data Summaries\sessionList.xlsx';
    [v, basepaths] = getSessionVars('basepaths', basepaths, 'varsFile', varsFile,...
        'varsName', varsName, 'pcond', ["tempflag"], 'ncond', [""],...
        'xlsname', xlsname);
    nfiles = length(basepaths);

    fr = catfields([v(:).fr], 'catdef', 'cell');

    for ifile = 1 : nfiles
        for iunit = 1 : 2
            unitIdx = v(ifile).units.clean(iunit, :);
            for istate = 1 : length(sstates)
                mfr(imouse, ifile, iunit, istate) = mean(fr.states.mfr{ifile}(unitIdx, sstates(istate)));
            end
        end
    end
end

squeeze(mfr(:, :, 2, 1))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MFR (in states) across mice per mouse
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CHECK SESSION_LIST.XLSX

sstates = [1, 4];
mname = {'lh123'; 'lh126'; 'lh129'; 'lh130'};
fileIdx = [2, 2, 4, 4];
binIdx = [1, 4, 1, 4];
clear mfr_states
for imouse = 1 : length(mname)
    varsFile = ["fr"; "fr_bins"; "datInfo"; "session"; "units"];
    varsName = ["fr"; "frBins"; "datInfo"; "session"; "units"];
    xlsname = 'D:\Google Drive\PhD\Slutsky\Data Summaries\sessionList.xlsx';
    [v, basepaths] = getSessionVars('mname', mname{imouse}, 'varsFile', varsFile,...
        'varsName', varsName, 'pcond', ["tempflag"], 'ncond', [""],...
        'xlsname', xlsname);
    nfiles = length(basepaths);


    for ibin = 1 : length(binIdx)
        mfrTemp_states = v(fileIdx(ibin)).frBins(binIdx(ibin)).states.mfr;
        mfrTemp = v(fileIdx(ibin)).frBins(binIdx(ibin)).mfr;
        for iunit = 1 : 2
            unitIdx = v(fileIdx(ibin)).units.clean(iunit, :);
            for istate = 1 : length(sstates)
                mfr_states(imouse, ibin, iunit, istate) = mean(mfrTemp_states(unitIdx, sstates(istate)), 'omitnan');
            end
            mfr(imouse, ibin, iunit) = mean(mfrTemp(unitIdx), 'omitnan');
        end
    end
end
m4_states = mfr_states;
m4 = mfr;

mname = {'lh122'};
fileIdx = [1, 2, 2, 3];
binIdx = [4, 1, 4, 1];
clear mfr_states mfr
for imouse = 1 : length(mname)
    varsFile = ["fr"; "fr_bins"; "datInfo"; "session"; "units"];
    varsName = ["fr"; "frBins"; "datInfo"; "session"; "units"];
    xlsname = 'D:\Google Drive\PhD\Slutsky\Data Summaries\sessionList.xlsx';
    [v, basepaths] = getSessionVars('mname', mname{imouse}, 'varsFile', varsFile,...
        'varsName', varsName, 'pcond', ["tempflag"], 'ncond', [""],...
        'xlsname', xlsname);
    nfiles = length(basepaths);


    for ibin = 1 : length(binIdx)
        mfrTemp_states = v(fileIdx(ibin)).frBins(binIdx(ibin)).states.mfr;
        mfrTemp = v(fileIdx(ibin)).frBins(binIdx(ibin)).mfr;
        for iunit = 1 : 2
            unitIdx = v(fileIdx(ibin)).units.clean(iunit, :);
            for istate = 1 : length(sstates)
                mfr_states(imouse, ibin, iunit, istate) = mean(mfrTemp_states(unitIdx, sstates(istate)), 'omitnan');
            end
            mfr(imouse, ibin, iunit) = mean(mfrTemp(unitIdx), 'omitnan');
        end
    end
end

iunit = 1;
istate = 1;
[squeeze(m4_states(:, :, iunit, istate)); squeeze(mfr_states(:, :, iunit, istate))]
[squeeze(m4(:, :, iunit)); squeeze(mfr(:, :, iunit))]


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get MFR per unit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mname = {'lh122'; 'lh123'; 'lh126'; 'lh129'; 'lh130'};

mfr = cell(length(sstates), 2, 3);
for imouse = 1 : length(mname)

    % reload data
    varsFile = ["fr"; "datInfo"; "session"; "units"];
    varsName = ["fr"; "datInfo"; "session"; "units"];
    xlsname = 'D:\Google Drive\PhD\Slutsky\Data Summaries\sessionList.xlsx';
    [v, basepaths] = getSessionVars('mname', mname{imouse}, 'varsFile', varsFile,...
        'varsName', varsName, 'pcond', ["tempflag"], 'ncond', [""],...
        'xlsname', xlsname);
    nfiles = length(basepaths);

    % organize in cell array
    for ifile = 1 : nfiles
        for iunit = 1 : 2
            for istate = 1 : length(sstates)
                unitIdx = v(ifile).units.clean(iunit, :);
                tmpMfr = v(ifile).fr.states.mfr(unitIdx, sstates(istate));
                mfr{istate, iunit, ifile} = [mfr{istate, iunit, ifile}; tmpMfr]
            end
        end
    end
end

% state ratio according to mfr percetiles
iunit = 2;
ifile = 1;
for ifile = 1 : 3
    plot_FRstates_sextiles('stateMfr', [mfr{:, iunit, ifile}]', 'units', [],...
        'ntiles', 4, 'saveFig', false)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% acsf single session
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

basepaths = [...
    "F:\Data\lh122\lh122_230109_095412";...
    "F:\Data\lh123\lh123_221228_102653";...
    "F:\Data\lh126\lh126_230117_102353";...
    "F:\Data\lh129\lh129_230214_093124";...
    "F:\Data\lh130\lh130_230405_094750"];

% reload data
varsFile = ["fr"; "fr_bins"; "datInfo"; "session"; "units"];
varsName = ["fr"; "frBins"; "datInfo"; "session"; "units"];
xlsname = 'D:\Google Drive\PhD\Slutsky\Data Summaries\sessionList.xlsx';
[v, ~] = getSessionVars('basepaths', basepaths, 'varsFile', varsFile,...
    'varsName', varsName, 'pcond', ["tempflag"], 'ncond', [""],...
    'xlsname', xlsname);
nfiles = length(basepaths);


% fr vs. time
iunit = 2;
[frMat, timeIdx] = alignFR2pnt('basepaths', basepaths, 'suFlag', true,...
    'dataType', 'strd', 'iunit', iunit, 'timeIdx', [0 0 0 0 0]);

% to prism
npnts = 13;
prismMat = [smooth(mean(frMat, 2), npnts),...
    smooth(std(frMat, [], 2) / sqrt(length(frMat)), npnts)];
nunits = sum(~isnan(frMat)');

% replace last hour with first hour
xval = [-6 : 1 / 60 : 18];
shiftVal = 90;  % [min]
[~, injIdx] = min(abs(xval - 0));
injIdx = find(injIdx);
prismMat(injIdx + shiftVal : end, :) = prismMat(injIdx : end - shiftVal, :)

% fr timebins per state
sstates = [1, 4];
mfr = cell(4, length(sstates), 2);
cnt = 1;
for ifile = 1 : nfiles
    for ibin = 1 : 4
        for iunit = 1 : 2
            unitIdx = v(ifile).units.clean(iunit, :);
            for istate = 1 : length(sstates)
                data = v(ifile).frBins(ibin).states.mfr(unitIdx, sstates(istate));
                mfr{ibin, iunit, istate} = [mfr{ibin, iunit, istate}; data];
            end
        end
    end
    cnt = cnt + 1;
end

iunit = 2;
istate = 2;
cell2nanmat(squeeze(mfr(:, iunit, istate)), 2)


% per mouse: MFR across units per state per unit
clear mfr
for ifile = 1 : nfiles
    for iunit = 1 : 2
        unitIdx = v(ifile).units.clean(iunit, :);
        for istate = 1 : 2
            mfr(1, ifile, iunit, istate) =...
                mean(v(ifile).frBins(1).states.mfr(unitIdx, sstates(istate)), 'omitnan')
             mfr(2, ifile, iunit, istate) =...
                mean(v(ifile).frBins(3).states.mfr(unitIdx, sstates(istate)), 'omitnan')
        end
    end
end

iunit = 2;
istate = 1;
squeeze(mfr(:, :, iunit, istate))'


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% acsf two sessions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -------------------------------------------------------------------------
% MFR vs time

clear basepaths
basepaths = ["F:\Data\lh107\lh107_220517_094900", "F:\Data\lh107\lh107_220518_091200";...
    "F:\Data\lh123\lh123_221228_102653", "F:\Data\lh123\lh123_221229_090102";...
    "F:\Data\lh126\lh126_230115_090453", "F:\Data\lh126\lh126_230117_102353";...
    "F:\Data\lh129\lh129_230214_093124", "F:\Data\lh129\lh129_230215_092653";...
    "F:\Data\lh130\lh130_230405_094750", "F:\Data\lh130\lh130_230406_091232"];
nmice = length(basepaths);
nfiles = numel(basepaths);
nsessions = size(basepaths, 2);

% fr vs. time
iunit = 2;
[frMat, timeIdx] = alignFR2pnt('basepaths', {basepaths{:, 1}}, 'suFlag', true,...
    'dataType', 'strd', 'iunit', iunit, 'timeIdx', [0 0 0 0 0]);
[frMat2, timeIdx] = alignFR2pnt('basepaths', {basepaths{:, 2}}, 'suFlag', true,...
    'dataType', 'strd', 'iunit', iunit, 'timeIdx', [0 0 0 0 0]);

npt = 91;
prismMfr = movmean([mean(frMat, 2, 'omitnan'); mean(frMat2, 2, 'omitnan')], npt);
prismSfr = movmean([std(frMat, [], 2, 'omitnan') / sqrt(size(frMat, 2));...
    std(frMat2, [], 2, 'omitnan') / sqrt(size(frMat, 2))], npt);
nunits = [ones(1, length(frMat)) * size(frMat, 2), ones(1, length(frMat2)) * size(frMat2, 2)];
xData = [1 : (length(frMat) + length(frMat2))] / 60 - 24
 
% -------------------------------------------------------------------------
% mfr in bins

% reload data
varsFile = ["fr"; "fr_bins"; "datInfo"; "session"; "units"];
varsName = ["fr"; "frBins"; "datInfo"; "session"; "units"];
xlsname = 'D:\Google Drive\PhD\Slutsky\Data Summaries\sessionList.xlsx';
[v, ~] = getSessionVars('basepaths', sort(basepaths(:)), 'varsFile', varsFile,...
    'varsName', varsName, 'pcond', ["tempflag"], 'ncond', [""],...
    'xlsname', xlsname);

% mfr per mouse between baseline and acute
iunit = 2;
clear mfr tmp
for ifile = 1 : nfiles

    unitIdx = v(ifile).units.clean(iunit, :);
    tmp(ifile) = mean(v(ifile).frBins(1).mfr(unitIdx));
    
end
mfr = [tmp(1 : 2 : end); tmp(2 : 2 : end)]';

% mfr per unit across bins
iunit = 1;
clear tmp
mfr = cell(1, 2);
for isession = 1 : nsessions
    for imouse = 1 : nmice

        unitIdx = v(imouse * isession).units.clean(iunit, :);
        tmp = [v(imouse * isession).frBins(:).mfr];
        mfr{isession} = [mfr{isession}; tmp(unitIdx, :)];

    end
end
prismMat = cell2nanmat(mfr);
mean(prismMat, 'omitnan')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cell classification
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mname = {'lh122'; 'lh123'; 'lh126'; 'lh129'; 'lh130'};
ifile = 1;

tp = []; brst = []; un = [];
for imouse = 1 : length(mname)

    % reload data
    varsFile = ["session"; "cell_metrics"; "units"];
    varsName = ["session"; "cm"; "units"];
    xlsname = 'D:\Google Drive\PhD\Slutsky\Data Summaries\sessionList.xlsx';
    [v, basepaths] = getSessionVars('mname', mname{imouse}, 'varsFile', varsFile,...
        'varsName', varsName, 'pcond', ["tempflag"], 'ncond', [""],...
        'xlsname', xlsname);
    nfiles = length(basepaths);
    
    for ifile = 1
        basepath = basepaths{ifile};
        cd(basepath)
        tp = [tp, v(ifile).cm.troughToPeak];
        brst = [brst, v(ifile).cm.st_royer];
        un = [un, v(ifile).units.clean];
    end

end

fh = figure;
iunit = 1;
prismMat = [tp(logical(un(iunit, :)))', brst(logical(un(iunit, :)))'];
scatter(prismMat(:, 1), prismMat(:, 2), 'k', 'filled')
hold on
iunit = 2;
prismMat = [tp(logical(un(iunit, :)))', brst(logical(un(iunit, :)))'];
scatter(prismMat(:, 1), prismMat(:, 2), 'r', 'filled')
set(gca, 'yscale', 'log')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% example waveform
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

basepath = 'F:\Data\lh129\lh129_230127_095804';
cd(basepath)
[~, basename] = fileparts(basepath);
cm = CellExplorer('basepath', pwd);
sessionFile = fullfile(basepath, [basename, '.session.mat']);
spkFile = fullfile(basepath, [basename, '.spikes.cellinfo.mat']);
load(sessionFile);
load(spkFile);


fs = 24;
rs = 35;
unitIdx = [rs, fs];

igrp = 3;
spk = loadNS('datatype', 'spk', 'session', session, 'grpid', igrp);
clu = loadNS('datatype', 'clu', 'session', session, 'grpid', igrp);
uclu = unique(clu);
spks2snip = 1000;

for iunit = 1 : 2
    
    % randomly select a subset of spikes
    cluId = spikes.cluID(spikes.UID == unitIdx(iunit));

    cluidx = find(clu == cluId);
    randidx = randperm(length(cluidx), min([length(cluidx), spks2snip]));
    cluidx = cluidx(randidx)
    
    wv{iunit} = spk(:, :, cluidx) * 0.195;
    
end

iunit = 2
fh= figure;
wvMean = mean(wv{iunit}, 3, 'omitnan');
wvSe = std(wv{iunit}, [], 3, 'omitnan') / sqrt(spks2snip);
nunits = ones(length(wvSe), 1)' * spks2snip;
plot_wv('wv', wvMean, 'wv_std', wvSe)
fs = 20000;
xval = ([-16 : 15] / fs) * 1000

cnt = 1;
for ich = 1 : 4
    prismMat(cnt : cnt + 2, :) = [wvMean(ich, :); wvSe(ich, :); nunits]
    cnt = cnt + 3;
end

% plots mean +- std waveform of all channels horizontal (concatenated) or
% vertical
%
% INPUT
%   wv          average waveform [nchans x nsamps] 
%   wv_std      std of waveform [nchams x nsamps]
%   clr         color of plot
%   orient      channels horizontal {horz} or vertical (vert)
%   sbar        plot scale bar {1} or not (0)
%   fs          sampling rate
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compare mfr (per unit) w/ msr (per tetrode)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


basepaths = mk801_chronic_sessions('bsl');

varsFile = ["session"; "fr"; "sr"; "units"];
varsName = ["session"; "fr"; "sr"; "units"];
xlsname = 'D:\Google Drive\PhD\Slutsky\Data Summaries\sessionList.xlsx';
[v, basepaths] = getSessionVars('basepaths', basepaths, 'varsFile', varsFile,...
    'varsName', varsName, 'pcond', ["tempflag"], 'ncond', [""],...
    'xlsname', xlsname);
nfiles = length(basepaths);

mfr = []; msr = [];
for ifile = 1 : nfiles
    
    unitIdx = v(ifile).units.clean(1, :);
    mfr = [mfr; v(ifile).fr.mfr(unitIdx)];
    msr = [msr; v(ifile).sr.mfr];

end

sqrt(var(mfr)) / mean(mfr)
sqrt(var(msr)) / mean(msr)



