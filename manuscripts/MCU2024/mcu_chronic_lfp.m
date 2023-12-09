

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

mname = {'lh96', 'lh107', 'lh122', 'lh142'};
mname = {'lh133'};
flgNormTime = true;
flgNormBand = false;

clear bands powdb
% bands(imouse, iband, istate, isession)
% powdb(imouse, istate, ifreq, isession)
for imouse = 1 : length(mname)
  
    [bands(imouse, :, :, :), powdb(imouse, :, :, :)] =...
        sessions_psd(mname{imouse}, 'flgNormBand', flgNormBand,...
        'flgAnalyze', false, 'flgNormTime', flgNormTime,...
        'flgEmg', true, 'idxBsl', [1], 'graphics', true);

end

% get mouse data (transpose s.t. sessions are rows). if flgNormBand, ignore
% first column (broadband)
istate = 2;
imouse = 1;
prismData = squeeze(bands(imouse, :, istate, :))';
prismData = squeeze(powdb(imouse, istate, :, :));

% transpose grp data
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

%%%% for each freq point, calc diff between bsl and bac effect, and plot as
%%%% function of freq point


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
% hynogram during baseline
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
varsFile = ["sleep_states"; "datInfo"; "session"];
varsName = ["ss"; "datInfo"; "session"];
xlsname = 'D:\Google Drive\PhD\Slutsky\Data Summaries\sessionList.xlsx';
[v, basepaths] = getSessionVars('basepaths', basepaths, 'varsFile', varsFile,...
    'varsName', varsName, 'pcond', ["tempflag"], 'ncond', [""],...
    'xlsname', xlsname);

% params
nfiles = length(basepaths);

% check states
for ifile = 1 : nfiles
    cd(basepaths{ifile})
    [~, basename] = fileparts(basepaths{ifile})
    sSig = load([basename, '.sleep_sig.mat']);
    AccuSleep_viewer(sSig, v(ifile).ss.labels, [])
end

% % add timebins to session
% for ifile = 1 : nfiles
%     cd(basepaths{ifile})
%     [timebins, timepnt] = metaInfo_timebins('reqPnt', 6 * 60 * 60, 'nbins', 4);
%     timebins / 60 / 60
% end


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
istate = 3;
prismData = prctDur(:, :, istate)';











