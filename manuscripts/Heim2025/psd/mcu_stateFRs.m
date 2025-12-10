
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FRs per state - calculate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load data
mname = 'mcu_wsh';
vars = ["sleep_states"; "session"; "units"; "spikes"; "psd"];
basepaths = mcu_basepaths(mname);
v = basepaths2vars('basepaths', basepaths, 'vars', vars);
nfiles = length(basepaths);

% params
sstates = [1, 4, 5];
clr = v(1).ss.info.colors(sstates);
snames = v(1).ss.info.names(sstates);
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
        'graphics', false, 'binsize', 60, 'saveVar', true, 'forceA', true,...
        'smet', 'none', 'winBL', [0, Inf], 'winCalc', [0, Inf],...
        'boutTimes', boutTimes);

    if graphics
        unitIdx = v(ifile).units.clean(iunit, :);
        plot_fr_boutTimes(fr, 'unitIdx', unitIdx, 'saveFig', false)
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FRs per state - load and organize
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% params
iunit = 2;
sstates = [1, 4, 5];

% grps
clear mname
mname{1} = 'wt_bsl';
mname{2} = 'mcu_bsl';
mname{1} = 'wt_wsh';
mname{2} = 'mcu_wsh';

clear boutGrp unitGrp
for igrp = 1 : length(mname)

    % load data
    vars = ["fr"; "units"];
    basepaths = mcu_basepaths(mname{igrp});
    v = basepaths2vars('basepaths', basepaths, 'vars', vars);
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
clear prismGrp
for igrp = 1 : 2;
    tmp = squeeze(unitMat(:, :, :, igrp));
    prismData = reshape(permute(tmp, [1, 3, 2]), [], 3);
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
igrp = 1;
tmp = squeeze(boutMat(:, :, :, igrp));
prismData = squeeze(mean(tmp, 1, 'omitnan'))

% Normalize MFR per mouse to AW
prismData = prismData ./ mean(prismData, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nunits
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% grps
clear mname
mname{1} = 'wt_bsl';
mname{2} = 'wt_wsh';
mname{3} = 'mcu_bsl';
mname{4} = 'mcu_wsh';
ngrp = length(mname)

nunits = nan(ngrp, 6, 2);
for igrp = 1 : ngrp

    % load data
    vars = ["units"];
    basepaths = mcu_basepaths(mname{igrp});
    v = basepaths2vars('basepaths', basepaths, 'vars', vars);
    nfiles = length(basepaths);

    % nunits
    for ifile = 1 : nfiles
        for iunit = 1 : 2
            unitIdx = v(ifile).units.clean(iunit, :);
            nunits(igrp, ifile, iunit) = sum(unitIdx);
        end
    end
end

iunit = 2;
sz = size(nunits);
prismData = squeeze(nunits([1 : 2], :, iunit));
prismData(:, sz(2) + 1 : sz(2) * 2) = squeeze(nunits([3, 4], :, iunit));

iunit = 2;
sum(nunits(:, :, iunit), 2, 'omitnan')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mfr per mouse irrespective of state
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% grps
clear mname
mname{1} = 'wt_bsl';
mname{2} = 'wt_wsh';
mname{3} = 'mcu_bsl';
mname{4} = 'mcu_wsh';
ngrp = length(mname)

mfr = nan(ngrp, 6, 2);
for igrp = 1 : ngrp

    % load data
    vars = ["fr"; "units"];
    basepaths = mcu_basepaths(mname{igrp});
    v = basepaths2vars('basepaths', basepaths, 'vars', vars);
    nfiles = length(basepaths);

    % firing rate
    for ifile = 1 : nfiles
        for iunit = 1 : 2
            unitIdx = v(ifile).units.clean(iunit, :);
            mfr(igrp, ifile, iunit) = mean(v(ifile).fr.mfr(unitIdx));
        end
    end
end

iunit = 1;
sz = size(mfr);
prismData = squeeze(mfr([1 : 2], :, iunit));
prismData(:, sz(2) + 1 : sz(2) * 2) = squeeze(mfr([3, 4], :, iunit));

iunit = 2;
sum(mfr(:, :, iunit), 2, 'omitnan')

