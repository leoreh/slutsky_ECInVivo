% wrapper for batch processing

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mname = 'lh96';
forceL = true;
forceA = true;

pcond = ["tempflag"];
ncond = [""];

% load vars from each session
varsFile = ["fr"; "sr"; "spikes"; "st_metrics"; "swv_metrics";...
    "cell_metrics"; "sleep_states"; "ripp.mat"; "datInfo"; "session"];
varsName = ["fr"; "sr"; "spikes"; "st"; "swv"; "cm"; "ss"; "ripp";...
    "datInfo"; "session"];
if ~exist('v', 'var') || forceL
    [v, basepaths] = getSessionVars('mname', mname, 'varsFile', varsFile,...
        'varsName', varsName, 'pcond', pcond, 'ncond', ncond);
end
nsessions = length(basepaths);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% analyze data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
templateCal = ss.info.calibrationData;

if forceA
    for isession = 1 : nsessions
        
        % file
        basepath = basepaths{isession};
        cd(basepath)
        [~, basename] = fileparts(basepath);
        
        % params
        session = v(isession).session;
        nchans = session.extracellular.nChannels;
        fs = session.extracellular.sr;
        spkgrp = session.extracellular.spikeGroups.channels;
        
        % add timebins to datInfo
        nbins = 4;
        reqPnt = 6 * 60 * 60;
        [timebins, timepnt] = metaInfo_timebins('reqPnt', reqPnt,...
            'nbins', nbins);
        winCalc = mat2cell(timebins, ones(nbins, 1), 2);
        
        % spk lfp
        frange = [0.5, 4; 5, 12; 50, 80];
        sl = spklfp_wrapper('basepath', basepath, 'winCalc', winCalc,...
            'ch', 9, 'frange', frange,...
            'graphics', true, 'saveVar', false);

        % spike timing metrics
        st = spktimesMetrics('winCalc', winCalc, 'forceA', true);


        % update units
        %         units = selectUnits('basepath', pwd, 'grp', [], 'saveVar', true,...
        %             'forceA', true, 'frBoundries', [0.1 Inf; 0.1 Inf],...
        %             'spikes', []);
           
        %         tbins_txt = {'0-3ZT', '3-6ZT', '6-9ZT', '9-12ZT',...
        %             '12-15ZT', '15-18ZT', '18-21ZT', '21-24ZT'};
        %         psdBins = psd_states_timebins('basepath', pwd,...
        %             'chEeg', [], 'forceA', true, 'graphics', true,...
        %             'timebins', chunks, 'saveVar', true,...
        %             'sstates', [1, 4, 5], 'tbins_txt', tbins_txt);
        
        %         psdBins = psd_timebins('basepath', pwd,...
        %             'chEeg', [34], 'forceA', true, 'graphics', true,...
        %             'timebins', chunks, 'saveVar', true);
        
        %         fr = firingRate(v(isession).spikes.times, 'basepath', basepath, 'graphics', true,...
        %             'binsize', 60, 'saveVar', true, 'smet', 'GK', 'winBL',...
        %             [0 timepoints], 'winCalc', [0, Inf], 'forceA', true);
        
        %         frBins(isession) = fr_timebins('basepath', pwd,...
        %             'forceA', false, 'graphics', true,...
        %             'timebins', chunks, 'saveVar', true);
        
    end
end

cell_metrics = CellExplorer('basepaths', basepaths);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% general params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

session = CE_sessionTemplate(pwd, 'viaGUI', false,...
    'force', true, 'saveVar', true);
nchans = session.extracellular.nChannels;
fs = session.extracellular.sr;
spkgrp = session.extracellular.spikeGroups.channels;

sessionIdx = 1 : nsessions;
stateidx = [1, 4, 5];
grp = [1 : 4];                  % which tetrodes to plot
unitClass = 'pyr';              % plot 'int', 'pyr', or 'all'
suFlag = 1;                     % plot only su or all units
frBoundries = [0 Inf];          % include only units with fr greater than

[nsub] = numSubplots(length(sessionIdx));
[cfg_colors, cfg_names, ~] = as_loadConfig([]);
setMatlabGraphics(false)

% arrange title names
for isession = 1 : nsessions
    sessionName{isession} = dirnames{isession}(length(mname) + 2 : end);
    basepath = char(fullfile(mousepath, dirnames{isession}));
    basepaths{isession} = fullfile(mousepath, dirnames{isession});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% states
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

