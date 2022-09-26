

% to prism

% load data
mname = 'lh96';

varsFile = ["fr"; "sr"; "spikes"; "st_metrics"; "swv_metrics";...
    "cell_metrics"; "sleep_states"; "ripp.mat"; "datInfo"; "session";...
    "psd_bins"; "units"];
varsName = ["fr"; "sr"; "spikes"; "st"; "swv"; "cm"; "ss"; "ripp";...
    "datInfo"; "session"; "psdBins"; "units"];
xlsname = 'D:\Google Drive\PhD\Slutsky\Data Summaries\sessionList.xlsx';
[v, basepaths] = getSessionVars('mname', mname, 'varsFile', varsFile,...
    'varsName', varsName, 'pcond', ["tempflag"], 'ncond', [""],...
    'xlsname', xlsname);
nfiles = length(basepaths);


% selections
ifile = 7;
sstates = [1, 4, 5];
iunit = 2;

% -------------------------------------------------------------------------
% single session
unitIdx = v(ifile).units.clean(iunit, :);

% fr vs. time
ydata = v(ifile).fr.strd(unitIdx, :)';
xdata = v(ifile).fr.tstamps / 60 / 60;


% -------------------------------------------------------------------------
% multiple sessions


% state ratio across sessions
fr = catfields([v(:).fr], 'catdef', 'cell');

clear sratio
for ifile = 1 : nfiles
    unitIdx = v(ifile).units.clean(iunit, :);
    sratio{ifile} = squeeze(fr.states.ratio{ifile}(4, 1, unitIdx));
    sgain{ifile} = squeeze(fr.states.gain{ifile}(4, unitIdx));
end
ydata = cell2nanmat(sratio, 2);
ydata = cell2nanmat(sgain, 2);


% how many units per session

clear clu4grp
for ifile = 1 : nfiles
    clu4grp(ifile, :) = histcounts(v(ifile).spikes.shankID);
end

ngrp = 4;
clear nunits
for ifile = 1 : nfiles
    for igrp = 1 : ngrp
        nunits(:, ifile, igrp) = sum(v(ifile).units.clean & v(ifile).spikes.shankID == igrp, 2);
    end
end
sum(squeeze(nunits(1, :, :)), 2)

