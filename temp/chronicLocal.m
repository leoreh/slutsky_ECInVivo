

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cfg = as_loadConfig;
sstates = [1, 4];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


mname = 'lh126';

varsFile = ["fr"; "sr"; "sleep_states";...
    "datInfo"; "session"; "units"; "psd"];
varsName = ["fr"; "sr"; "ss"; "datInfo"; "session";...
    "units"; "psd"];
xlsname = 'D:\Google Drive\PhD\Slutsky\Data Summaries\sessionList.xlsx';
[v, basepaths] = getSessionVars('mname', mname, 'varsFile', varsFile,...
    'varsName', varsName, 'pcond', ["tempflag"], 'ncond', [""],...
    'xlsname', xlsname);
nfiles = length(basepaths);



% get psd per session for one mouse
bands = sessions_psd(mname, 'flgNormBand', true, 'flgAnalyze', false,...
    'flgNormTime', true, 'flgEmg', true, 'idxBsl', [1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% organize data - spiking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fr = [v(:).fr];
states = catfields([fr(:).states], 'catdef', 'cell', 'force', false);



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


