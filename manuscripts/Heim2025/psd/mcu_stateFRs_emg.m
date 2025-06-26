
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FRs per state - calculate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load data
mname = 'wt_bsl';
vars = ["sleep_states"; "session"; "units"; "spikes"; "psdEmg"];
[basepaths, v] = mcu_sessions(mname, vars);
nfiles = length(basepaths);

% params
sstates = [1, 2];
clr = v(1).psd.info.clr;
snames = v(1).psd.info.snames;
graphics = true;
yLimit_fr = [0, 3];
iunit = 1;

% firing rate
for ifile = 1 : nfiles

    basepath = basepaths{ifile};
    [~, basename] = fileparts(basepath);
    cd(basepath)

    % get state bouts from psd
    boutTimes = v(ifile).psd.bouts.boutTimes;
    %     boutTimes = v(ifile).ss.boutTimes(sstates);

    % get mfr per bout
    fr = calc_fr(v(ifile).spikes.times, 'basepath', basepath,...
        'graphics', false, 'binsize', 60, 'saveVar', 'frEmg', 'forceA', true,...
        'smet', 'none', 'winBL', [0, Inf], 'winCalc', [0, Inf],...
        'boutTimes', boutTimes);

    if graphics
        unitIdx = v(ifile).units.clean(iunit, :);
        plot_fr_boutTimes(fr(ifile), 'unitIdx', unitIdx,...
            'sstates', sstates, 'saveFig', false)
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FRs per state - load and organize
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% params
iunit = 2;
sstates = [1, 2];

% grps
clear mname
mname{1} = 'wt_bsl';
mname{2} = 'mcu_bsl';
% mname{1} = 'wt_wsh';
% mname{2} = 'mcu_wsh';

clear boutGrp unitGrp
for igrp = 1 : length(mname)

    % load data
    vars = ["frEmg"; "units"];
    [basepaths, v] = mcu_sessions(mname{igrp}, vars);
    nfiles = length(basepaths);

    % firing rate
    clear unitCell boutCell
    for ifile = 1 : nfiles

        % get mfr across units per state bout
        unitIdx = v(ifile).units.clean(iunit, :);
        tmp = cellfun(@(x) mean(x(unitIdx, :), 1)', v(ifile).fr.states.fr, 'uni', false);
        tmp = cell2padmat(tmp, 2);
        boutCell{ifile} = tmp;

        % get mfr across state bouts per unit
        unitCell{ifile} = v(ifile).fr.states.mfr(unitIdx, :);        
        
    end
    boutGrp{igrp} = cell2padmat(boutCell, 3);
    unitGrp{igrp} = cell2padmat(unitCell, 3);

end
boutMat = cell2padmat(boutGrp, 4);
unitMat = cell2padmat(unitGrp, 4);

% FRs per unit across bouts, all units
clear prismGrp prismData
for igrp = 1 : 2;
    tmp = squeeze(unitMat(:, :, :, igrp));
    prismData = reshape(permute(tmp, [1, 3, 2]), [], length(sstates));
    nunits = sum(~isnan(prismData));
    prismGrp{igrp} = [mean(prismData, 'omitnan')',...
        std(prismData, 'omitnan')' / sqrt(length(prismData)),...
        nunits'];
end
prismData = cat(2, prismGrp{:});

% MFR per bouts across units, all bouts
clear prismGrp
for igrp = 1 : 2;
    tmp = squeeze(boutMat(:, :, :, igrp));
    prismData = reshape(permute(tmp, [1, 3, 2]), [], 3);
    nunits = sum(~isnan(prismData));
    prismGrp{igrp} = [mean(prismData, 'omitnan')',...
        std(prismData, 'omitnan')' / sqrt(length(prismData)),...
        nunits'];
end
prismData = cat(2, prismGrp{:});

% MFR per mouse
igrp = 2;
tmp = squeeze(boutMat(:, :, :, igrp));
prismData = squeeze(mean(tmp, 1, 'omitnan'))

% Normalize MFR per mouse to overall MFR
prismData = prismData ./ mean(prismData, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% state (emg) dependent MFR across the experiment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for each session (and mouse), calculate EMG states
% calculate fr in timebins

% gen parans
saveVar = false;
graphics = true;

% state params
minDur = 10;
interDur = 3;
ftarget = [0.5 : 0.5 : 100];

% files
mname = 'lh107';
basepaths = [mcu_sessions(mname)];
nfiles = length(basepaths);

for ifile = 1 : nfiles

    basepath = basepaths{ifile};
    [~, basename] = fileparts(basepath);
    cd(basepath)


    bouts = as_bouts('minDur', minDur, 'interDur',...
        interDur, 'flgEmg', true, 'graphics', graphics, 'nbins', 2);
    boutTimes = bouts.times;
    
    % calc psd according to emg state separation
    psd = psd_states('basepath', basepath, 'sstates', [1, 2],...
        'sig', [], 'fs', fs, 'saveVar', saveVar,...
        'graphics', true, 'forceA', true, 'ftarget', ftarget,...
        'emgThr', [], 'flgEmg', true, 'boutTimes', boutTimes);

end
