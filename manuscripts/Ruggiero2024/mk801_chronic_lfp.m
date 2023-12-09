

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cfg = as_loadConfig;
sstates = [1, 4];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


mname = 'lh134';

varsFile = ["fr"; "sr"; "sleep_states";...
    "datInfo"; "session"; "units"; "psd"];
varsName = ["fr"; "sr"; "ss"; "datInfo"; "session";...
    "units"; "psd"];
xlsname = 'D:\Google Drive\PhD\Slutsky\Data Summaries\sessionList.xlsx';
[v, basepaths] = getSessionVars('mname', mname, 'varsFile', varsFile,...
    'varsName', varsName, 'pcond', ["tempflag"], 'ncond', [""],...
    'xlsname', xlsname);
nfiles = length(basepaths);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% average bands and psd from multiple mice 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mname = {'lh96', 'lh107', 'lh122'};
% mname = {'lh132', 'lh133', 'lh134'};
flgNormTime = true;
flgNormBand = true;

clear bands powdb
for imouse = 1 : length(mname)
  
    [bands(imouse, :, :, :), powdb(imouse, :, :, :)] =...
        sessions_psd(mname{imouse}, 'flgNormBand', flgNormBand,...
        'flgAnalyze', false, 'flgNormTime', flgNormTime,...
        'flgEmg', true, 'idxBsl', [1], 'graphics', true);

end

istate = 1;
% transpose
prismData = squeeze(mean(bands(:, :, istate, :), 1, 'omitnan'));
prismData = squeeze(mean(powdb(:, istate, :, :), 1, 'omitnan'));

% reshape 3d to 2d for grouped graph in prism
x = squeeze(powdb(:, istate, :, :));
prismData = [];
cnt = 1;
for ifile = 1 : size(x, 3)
    for imouse = 1 : length(mname)
        prismData(:, cnt) = squeeze(x(imouse, :, ifile));
        cnt = cnt + 1;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% organize data - psd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

psd = catfields([v(:).psd], 'catdef', 'cell', 'force', false);

% ORGANIZE: statePsd is a matrix of frequency (rows) x session (column)
% depicting the psd for the selected state
clear statePsd 
istate = 4;    
for ifile = 1 : nfiles

    statePsd(:, ifile) = squeeze(psd.psd{ifile}(1, istate, :));

end
statePsd = statePsd ./ sum(statePsd);
freq = psd.info.faxis{1};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% organize data - ripples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iunit = 1;
clear rippGain
for ifile = 1 : nfiles

    unitIdx = v(ifile).units.clean(iunit, :);
    rippMfr = v(ifile).ripp.spks.su.rippMap(unitIdx, :, :);
    rippMfr = squeeze(mean(mean(rippMfr, 2), 3));
    randMfr = v(ifile).ripp.spks.su.randMap(unitIdx, :, :);
    randMfr = squeeze(mean(mean(randMfr, 2), 3));    
    rippGain{ifile} = (rippMfr - randMfr) ./ (rippMfr + randMfr);

end
rippGain = cell2nanmat(rippGain, 2);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% hynogram during baseline
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% mcu baseline basepaths
queryStr = 'mk801';
basepaths = mk801_chronic_sessions(queryStr)
nfiles = length(basepaths);

% load data
varsFile = ["sleep_states"; "datInfo"; "session"];
varsName = ["ss"; "datInfo"; "session"];
xlsname = 'D:\Google Drive\PhD\Slutsky\Data Summaries\sessionList.xlsx';
[v, basepaths] = getSessionVars('basepaths', basepaths, 'varsFile', varsFile,...
    'varsName', varsName, 'pcond', ["tempflag"], 'ncond', [""],...
    'xlsname', xlsname);

% check states
for ifile = 1 : nfiles
    cd(basepaths{ifile})
    [~, basename] = fileparts(basepaths{ifile})
    sSig = load([basename, '.sleep_sig.mat']);
    AccuSleep_viewer(sSig, v(ifile).ss.labels, [])
end

% add timebins to session
for ifile = 1 : nfiles
    cd(basepaths{ifile})
    [timebins, timepnt] = metaInfo_timebins('reqPnt', 6 * 60 * 60, 'nbins', 4);
    timebins / 60 / 60
end


% plot hynogram
fh = figure;
for ifile = 1 : nfiles
    
    sb = subplot(nfiles, 1, ifile);
    plot_hypnogram('stateEpochs', v(ifile).ss.stateEpochs, 'axh', sb)

end

% plot state duration in timebins
sstates = [1, 4, 5];
clear prctDur
for ifile = 1 : nfiles
    cd(basepaths{ifile})
    [totDur, prctDur(ifile, :, :), epLen] = as_plotZT('nwin', 4,...
        'sstates', sstates, 'ss', v(ifile).ss,...
        'timebins', v(ifile).session.general.timebins, 'graphics', false);
end

% reorganize for prism
istate = 2;
prismData = prctDur(:, :, istate);



