

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% analysis per mouse
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mname = 'lh107';

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
    fprintf('working on session %d of %d, %s\n',...
        ifile, nfiles, basename)

    % add timebins to datInfo
    [timebins, timepnt] = metaInfo_timebins('reqPnt', 6 * 60 * 60, 'nbins', 4);
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
    fprintf('working on session %d of %d, %s\n',...
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

mname = {'lh96'; 'lh107'; 'lh122'; 'lh142'};
mname = {'lh132'; 'lh133'; 'lh134'; 'lh136'; 'lh140'};

for imouse = 1 : length(mname)

    mfrcat = cell(3, 2);    % {states, unitType}; mfrcat{3, :} = mfr irrespective of states
    frMed = [];      % take units with mfr > med (pos), < med (neg), or all []
    sstates = [1, 4];

    % reload data
    varsFile = ["fr"; "fr_bins"; "datInfo"; "session"; "units"];
    varsName = ["fr"; "frBins"; "datInfo"; "session"; "units"];
    xlsname = 'D:\Google Drive\PhD\Slutsky\Data Summaries\sessionList.xlsx';
    [v, basepaths] = getSessionVars('mname', mname{imouse}, 'varsFile', varsFile,...
        'varsName', varsName, 'pcond', ["tempflag"], 'ncond', [""],...
        'xlsname', xlsname);
    nfiles = length(basepaths);

    % organize in cell array
    cnt = 1;
    clear mfr stateRat stateGain

    for ifile = 1 : nfiles
        for ibin = 1 : 4
            for iunit = 1 : 2

                unitIdx = v(ifile).units.clean(iunit, :);
                unitMfr = v(ifile).frBins(ibin).mfr(unitIdx)';
                allMfr = v(ifile).frBins(ibin).mfr';
                if isempty(frMed)
                    unitMfrIdx = ones(1, length(allMfr));
                elseif frMed > 0
                    unitMfrIdx = allMfr > median(unitMfr);
                elseif frMed < 0
                    unitMfrIdx = allMfr < median(unitMfr)';
                end
                unitIdx = unitIdx & unitMfrIdx;

                for istate = 1 : length(sstates)

                    mfr{cnt, istate, iunit} = v(ifile).frBins(ibin).states.mfr(unitIdx, sstates(istate));

                end
                mfr{cnt, 3, iunit} = v(ifile).frBins(ibin).mfr(unitIdx);
                stateRat{cnt, iunit} = squeeze(v(ifile).frBins(ibin).states.ratio(1, 4, unitIdx));
                stateGain{cnt, iunit} = squeeze(v(ifile).frBins(ibin).states.gain(4, unitIdx));
            end
            cnt = cnt + 1;
        end
    end

    % reorganize for prism
    for istate = 1 : 3
        for iunit = 1 : 2
            data = cell2nanmat(squeeze(mfr(:, istate, iunit)), 2);
            mfrcat{istate, iunit} = [mfrcat{istate, iunit}; data];
        end
    end
    iunit = 1;
    cell2nanmat(squeeze(stateRat(:, iunit)), 2);
    cell2nanmat(squeeze(stateGain(:, iunit)), 2);

    mfrcat{3, 2}

    % save to mousepath
    mpath = fullfile('F:', 'Data', mname(imouse));
    fname = fullfile(mpath{1}, [mname{imouse}, '.mfrcat.mat']);
    save(fname, 'mfrcat')

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get MFR per unit across mice (by experimental groups)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


mname = {'lh96'; 'lh107'; 'lh122'; 'lh142'};
mname = {'lh132'; 'lh133'; 'lh134'; 'lh136'; 'lh140'};
nmice = length(mname);

clear mfr
for imouse = 1 : nmice
    
    mpath = fullfile('F:', 'Data', mname(imouse));
    fname = fullfile(mpath{1}, [mname{imouse}, '.mfrcat.mat']);
    load(fname)
    mfr(imouse, :, :) = mfrcat;

end

istate = 3;         % mfr irrespective of states
iunit = 2;

clear mfrPrism
data = mfr(:, istate, iunit);
mfrPrism = data{1};
cnt = size(mfrPrism, 1);
for imouse = 2 : nmice
    nunits = size(data{imouse}, 1)
    mfrPrism(cnt + 1 : cnt + nunits, :) = data{imouse}
    cnt = cnt + nunits;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get MFR per mouse (by experimental groups)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% WT
mname = {'lh96'; 'lh107'; 'lh122'; 'lh142'};
% choose bins per mouse
bins(1, :) = [4, 6, 21, 22];
bins(2, :) = [4, 6, 20, 25];
bins(3, :) = [4, 7, 20, 24];
bins(4, :) = [2, 9, 20, 24];

mname = {'lh132'; 'lh133'; 'lh134'; 'lh136'; 'lh140'};
% choose bins per mouse
bins(1, :) = [4, 7, 20, 24];
bins(2, :) = [4, 9, 20, 23];
bins(3, :) = [4, 7, 20, 24];
bins(4, :) = [4, 9, 20, 24];
bins(5, :) = [4, 6, 20, 23];

nmice = length(mname);

istate = 3;         % mfr irrespective of states
iunit = 2;

clear mfr
for imouse = 1 : nmice
    
    mpath = fullfile('F:', 'Data', mname(imouse));
    fname = fullfile(mpath{1}, [mname{imouse}, '.mfrcat.mat']);
    load(fname)

    mfr(imouse, :) = mean(mfrcat{istate, iunit}(:, bins(imouse, :)), 'omitnan');

end

% normalize to baseline
mfr ./ mfr(:, 1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MFR per unit in states during baseline
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% mcu baseline basepaths
basepaths = {...
    'F:\Data\lh132\lh132_230413_094013',...
    'F:\Data\lh133\lh133_230413_094013',...
    'F:\Data\lh134\lh134_230504_091744',...
    'F:\Data\lh136\lh136_230519_090043',...
    'F:\Data\lh140\lh140_230619_090023'};

% wt baseline basepaths
basepaths = {...
    'F:\Data\lh96\lh96_220120_090157',...
    'F:\Data\lh107\lh107_220518_091200',...
    'F:\Data\lh122\lh122_221223_092656',...
    'F:\Data\lh142\lh142_231005_091832'};


% load data
varsFile = ["fr"; "sleep_states"; "datInfo"; "session"; "units";...
    "st_metrics"; "st_brst"; "spikes"];
varsName = ["fr"; "ss"; "datInfo"; "session"; "units";...
    "st"; "brst"; "spikes"];
xlsname = 'D:\Google Drive\PhD\Slutsky\Data Summaries\sessionList.xlsx';
[v, basepaths] = getSessionVars('basepaths', basepaths, 'varsFile', varsFile,...
    'varsName', varsName, 'pcond', ["tempflag"], 'ncond', [""],...
    'xlsname', xlsname);
nfiles = length(basepaths);

sstates = [1, 4, 5];
iunit = 1;

% mfr
data = [];
pmfr = nan(nfiles, length(sstates));
for ifile = 1 : nfiles

    unitIdx = v(ifile).units.clean(iunit, :);
    
    data = [data; v(ifile).fr.states.mfr(unitIdx, sstates)];
    pmfr(ifile, :) = mean(v(ifile).fr.states.mfr(unitIdx, sstates));

end

% burstiness from st
data = [];
dataType = 'royer';
for ifile = 1 : nfiles

    unitIdx = v(ifile).units.clean(iunit, :);
    
    % mfr
    data = [data, v(ifile).st.(dataType)(unitIdx)];
    
end

% burstiness from mea_brst
data = [];
dataType = 'rate';
for ifile = 1 : nfiles

    unitIdx = v(ifile).units.clean(iunit, :);
    
    % mfr
    data = [data, v(ifile).brst.(dataType)(unitIdx)];
    
end



%%% analyze
for ifile = 1 : nfiles
    cd(basepaths{ifile})
    brst = spktimes_meaBrst(v(ifile).spikes.times, 'binsize', [], 'isiThr', 0.05,...
        'minSpks', 2, 'saveVar', true, 'force', true, 'bins', [0 Inf]);
end



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


