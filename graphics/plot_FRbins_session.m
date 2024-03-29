
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% selection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

basepath = pwd;
[~, basename] = fileparts(basepath);
cd(basepath)
guessDateTime(basename)

saveFig = false;
grp = [2 : 4];                  % which tetrodes to plot
suFlag = true;                  % plot only su or all units
stateIdx = [1, 4];
% include only units with fr greater / lower than. 1st row RS 2nd row FS
frBoundries = [0.2 Inf; 0.2 Inf];  
cellType = ["RS"; "FS"];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[cfg_colors, cfg_names, ~] = as_loadConfig([]);
[varArray, ~, mousepath] = getSessionVars('dirnames', string(basename));
assignVars(varArray, 1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% select time bins
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% time bins are created in binsize steps from a specific point of interest
% or from the center of the recording <csec(end) / 2> to the start and end
% of the recording. thus, binsize = Inf will yield two time bins. 

fs = session.extracellular.sr;
binsize = 3 * 60 * 60;      % hr in [s]
binsize = Inf;              % will separate recording into two parts
csec = floor(cumsum(datInfo.nsamps) / fs);

% point of interest
expPoint = csec(end) / 2;
expPoint = csec(1);

% time bins
timebins = sort([expPoint : - binsize : 1,...
    expPoint + binsize : binsize : csec(end)]);
timebins = [[1, timebins + 1]; [timebins, csec(end)]]';

% validate bin length above thr
binsize_len = timebins(:, 2) - timebins(:, 1);
binsize_thr = binsize / 3;
if binsize < Inf
    timebins(binsize_len < binsize_thr, :) = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arrange data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[mfrCell, gainCell] = org_mfrCell('spikes', spikes, 'cm', cm, 'fr', fr,...
    'timebins', timebins, 'dataType', 'su', 'grp', grp, 'suFlag', suFlag,...
    'frBoundries', frBoundries, 'stateIdx', stateIdx);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% mfr
fh = figure;
count = 1;
for itype = 1 : 2
    for istate = 1 : length(stateIdx)
        subplot(2, length(stateIdx), count)
        mfrMat = mfrCell{itype, istate};
        plot(mfrMat')
        
        if length(timebins) == 2
            pval  = signrank(mfrMat(:, 1), mfrMat(:, 2));
        else
            pval = friedman(mfrMat, 1, 'off');
        end
        
        xlim([1 - 0.5, length(timebins) + 0.5])
        xticks(1 : length(timebins))
        set(gca, 'box', 'off', 'TickLength', [0 0])
        ylabel('MFR [Hz]')
        title(sprintf('%s in %s; p = %.2f',...
            cellType{itype}, cfg_names{stateIdx(istate)}, pval))
        count = count + 1;
    end
end

% gain
fh = figure;
for itype = 1 : 2
        subplot(2, 1, itype)
        gainMat = gainCell{itype};
        plot(gainMat')
        
        if length(timebins) == 2
            pval  = signrank(gainMat(:, 1), gainMat(:, 2));
        else
            pval = friedman(gainMat, 1, 'off');
        end
        
        xlim([1 - 0.5, length(timebins) + 0.5])
        xticks(1 : length(timebins))
        set(gca, 'box', 'off', 'TickLength', [0 0])
        ylabel('Gain Factor [%]')
        title(sprintf('%s; p = %.2f',...
            cellType{itype}, pval))
end


