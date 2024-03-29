

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


% fr in time bins bins
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
% get MFR in time bins per unit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% sessionList: all five sessions should be considered

mname = {'lh122'; 'lh123'; 'lh126'; 'lh129'; 'lh130'};

sstates = [1, 4];
mfrcat = cell(length(sstates), 2);
gaincat = cell(1, 2);
frMed = [-1];         % take units with mfr > med (pos), < med (neg), or all []
gainVar = 'ratio';   % can be 'gain' or 'ratio'
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
    clear mfr stateGain
    for ifile = 1 : nfiles
        for ibin = 1 : 4
            for iunit = 1 : 2

                unitIdx = v(ifile).units.clean(iunit, :);
                unitMfr = v(ifile).frBins(ibin).mfr(unitIdx)';
                allMfr = v(ifile).frBins(ibin).mfr';
                if isempty(frMed)
                    unitMfrIdx = ones(1, length(allMfr));
                elseif frMed > 0
                    unitMfrIdx = allMfr > prctile(unitMfr, 75);
                elseif frMed < 0
                    unitMfrIdx = allMfr < prctile(unitMfr, 25)';
                end
                unitIdx = unitIdx & unitMfrIdx;

                for istate = 1 : length(sstates)
                                       
                    mfr{cnt, istate, iunit} = v(ifile).frBins(ibin).states.mfr(unitIdx, sstates(istate));

                end
                
                % get state ratio / gain factor
                tmp = v(ifile).frBins(ibin).states.(gainVar);
                if ndims(tmp) == 3
                    stateGain{cnt, iunit} = squeeze(tmp(4, 1, unitIdx));
                elseif ndims(tmp) == 2
                    stateGain{cnt, iunit} = squeeze(tmp(4, unitIdx));
                end
            end
            cnt = cnt + 1;
        end
    end

    % reorganize for prism
    for istate = 1 : length(sstates)
        for iunit = 1 : 2
            data = cell2nanmat(squeeze(mfr(:, istate, iunit)), 2);
            mfrcat{istate, iunit} = [mfrcat{istate, iunit};...
                data];
        end
    end
    for iunit = 1 : 2
        data = cell2nanmat(squeeze(stateGain(:, iunit)), 2);
        gaincat{iunit} = [gaincat{iunit}; data];
    end

end

% to prism: mfr in bins 
iunit = 1;
mfrcat{2, iunit}

% to prism: scatter data of mfr wake vs mfr nrem vs gain factor
iunit = 1;
bins = -30 : 6 : 84;
clear data
for ibin = 1 : size(gaincat{1}, 2)
    data = [mfrcat{1, iunit}(:, ibin),...
        mfrcat{2, iunit}(:, ibin),...
        gaincat{iunit}(:, ibin)];

    % calc fraction of units with gain factor >0
    tmp = gaincat{iunit}(:, ibin);
    frct(ibin) = sum(tmp > 0) / length(tmp);
end

% correct lh122_230112_084703
% nunits = 27;
% for ibin = 1 : 4
%     frBins(ibin).mfr = nan(nunits, 1);
%     frBins(ibin).states.mfr = nan(nunits, 6);
%     frBins(ibin).states.gain = nan(6, nunits);
%     frBins(ibin).states.ratio = nan(6, 6, nunits);
% end
% basepath = pwd;
% [~, basename] = fileparts(basepath)
% filename = [basename, '.fr_bins.mat'];
% save(filename, 'frBins')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MFR in states across mice per unit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CHECK SESSION_LIST.XLSX

sstates = [1 : 6];
mname = {'lh122'; 'lh123'; 'lh126'; 'lh129'; 'lh130'};
clear mfr
for imouse = 1 : length(mname)
    varsFile = ["fr"; "datInfo"; "session"; "units"];
    varsName = ["fr"; "datInfo"; "session"; "units"];
    xlsname = 'D:\Google Drive\PhD\Slutsky\Data Summaries\sessionList.xlsx';
    [v, basepaths] = getSessionVars('mname', mname{imouse}, 'varsFile', varsFile,...
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
% MFR vs. time across mice (baclofen)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mname = {'lh96'; 'lh107'; 'lh122'};

fr = cell(length(mname), 2);
for imouse = 1 : length(mname)
    
    for iunit = 1 : 2
    
        % concatenate var from different sessions
        [expData, xData] = sessions_catVarTime('mname', mname{imouse},...
            'dataPreset', 'fr', 'graphics', false, 'dataAlt', iunit,...
            'basepaths', {}, 'xTicksBinsize', 6, 'markRecTrans', true);

        % remove nan from columns (units)
        expData(:, all(isnan(expData), 1)) = [];

        % put in cell
        fr{imouse, iunit} = expData;

    end
end

% smooth each cell and than calc mean and std and smooth across population

iunit = 1;
npt = 91;
prismMat = cell2nanmat(fr(:, iunit));
prismMat = movmean(prismMat, npt, 1, 'omitnan');
nunits = sum(~isnan(prismMat'))';
xData = [1 / 60 : 1 / 60 : length(nunits) / 60] - 30;

mfr = movmean(mean(prismMat, 2, 'omitnan'), npt);
sfr = movmean(std(prismMat, [], 2, 'omitnan') ./ sqrt(nunits), npt);
sfr(sfr == Inf) = nan;

fh = figure
plot(xData, mfr)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MFR per unit in time bins (baclofen / mcu-ko)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mname = {'lh96'; 'lh107'; 'lh122'};

iunit = 1;
mfr = cell(6, 1);
for imouse = 1 : length(mname)

    % reload data
    varsFile = ["fr"; "datInfo"; "session"; "units"; "fr_bins"];
    varsName = ["fr"; "datInfo"; "session"; "units"; "frBins"];
    xlsname = 'D:\Google Drive\PhD\Slutsky\Data Summaries\sessionList.xlsx';
    [v, basepaths] = getSessionVars('mname', mname{imouse}, 'varsFile', varsFile,...
        'varsName', varsName, 'pcond', ["tempflag"], 'ncond', [""],...
        'xlsname', xlsname);
    nfiles = length(basepaths);

    for ifile = 1 : nfiles

            unitIdx = v(ifile).units.clean(iunit, :);
            tmpStruct = catfields([v(ifile).frBins], 'catdef', 'addim');
            mfr{ifile} = [mfr{ifile}; tmpStruct.mfr(unitIdx, :)];
    end
end

prismMat = cell2nanmat(mfr);

fh = figure;
plot_boxMean('dataMat', prismMat)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MFR per unit in time bins (baclofen mcu-ko)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mname = {'lh132'; 'lh133'; 'lh134'};
npt = 131;

sr = cell(length(mname), 1);
for imouse = 1 : length(mname)
    
    
        % concatenate var from different sessions
        [expData, xData] = sessions_catVarTime('mname', mname{imouse},...
            'dataPreset', 'sr', 'graphics', false, 'dataAlt', 1,...
            'basepaths', {}, 'xTicksBinsize', 6, 'markRecTrans', true);

        % put in cell
        sr{imouse} = expData;
        srMean{imouse} = mean(movmean(expData, npt, 1, 'omitnan'), 2);
        
end

prismMat = cell2nanmat(sr, 1);
prismMat = movmean(prismMat, npt, 1, 'omitnan');


% recalculate spike rate in bins
sr = [];
for imouse = 1 : length(mname)
 
    % reload data
    varsFile = ["spktimes"; "session"];
    varsName = ["spktimes"; "session"];
    xlsname = 'D:\Google Drive\PhD\Slutsky\Data Summaries\sessionList.xlsx';
    [v, basepaths] = getSessionVars('mname', mname{imouse}, 'varsFile', varsFile,...
        'varsName', varsName, 'pcond', ["tempflag"], 'ncond', [""],...
        'xlsname', xlsname);
    nfiles = length(basepaths);
    
    srTmp2 = []
    for ifile = 1 : nfiles
        
        binsize = 60 * 60 * 6 * 20000;
        nsamps = v(ifile).session.general.nsamps;
        binedges = [0 : nsamps / 4 : nsamps];
        binedges(end) = nsamps;
        if length(binedges) ~= 5
            warning('stop')
        end
        
        srTmp = [];
        for itet = 1 : 4
            srTmp(itet, :) = histcounts(v(ifile).spktimes{itet}, binedges,...
                'Normalization', 'countdensity')' * 20000;
        end
        
        srTmp2 = [srTmp2, srTmp];

    end

    sr = [sr; srTmp2];

        
end

% mean per mouse
cnt = 0;
clear srMean
for imouse = 1 : length(mname)
    tetIdx = [1 : 4] + 4 * cnt;
    srMean(imouse, :) = mean(sr([tetIdx], [1, 9, 21, 28]));
    cnt = cnt + 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MFR per mouse during timepoints (baclofen)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mname = {'lh96'; 'lh107'; 'lh122'};

iunit = 1;
mfr = cell(length(mname), 1);
for imouse = 1 : length(mname)

    % reload data
    varsFile = ["fr"; "datInfo"; "session"; "units"; "fr_bins"];
    varsName = ["fr"; "datInfo"; "session"; "units"; "frBins"];
    xlsname = 'D:\Google Drive\PhD\Slutsky\Data Summaries\sessionList.xlsx';
    [v, basepaths] = getSessionVars('mname', mname{imouse}, 'varsFile', varsFile,...
        'varsName', varsName, 'pcond', ["tempflag"], 'ncond', [""],...
        'xlsname', xlsname);
    nfiles = length(basepaths);

    for ifile = 1 : nfiles

            unitIdx = v(ifile).units.clean(iunit, :);
            tmpStruct = catfields([v(ifile).frBins], 'catdef', 'addim');
            mfr{imouse} = [mfr{imouse}; mean(tmpStruct.mfr(unitIdx, :), 1, 'omitnan')'];

    end
end

prismMat = cell2nanmat(mfr, 2);
prismMat = prismMat([1, 7, 22, 24], :)';


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



