function ripp = rippleSpks(ripp, varargin)

% gets spiking activity in relation to ripples. mainly creates maps: per
% ripple for multiunit activity, and per unit for single unit firing rate.
% complemantary to getRipples.m - requires the struct ripp. can also
% implement the analysis from bz_getRipSpikes
%
% INPUT:
%   basepath            path to recording {pwd}
%   graphics            logical. plot graphics {true} or not (false)
%   saveVar             logical. save variables (update ripp)
%   fullAnalysisFlag    logical. implement bz_getRipSpikes. really time
%                       consuming
%
% OUTPUT:
%   ripp            struct
%
% DEPENDENCIES:
%   bz_getRipSPikes (buzcode)
%   plot_rippleSpks
%
% 12 jan 23 LH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addOptional(p, 'basepath', pwd);
addOptional(p, 'graphics', true, @islogical);
addOptional(p, 'saveVar', true, @islogical);
addOptional(p, 'fullAnalysisFlag', false, @islogical);

parse(p, varargin{:})
basepath            = p.Results.basepath;
graphics            = p.Results.graphics;
saveVar             = p.Results.saveVar;
fullAnalysisFlag    = p.Results.fullAnalysisFlag;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% files
cd(basepath)
[~, basename] = fileparts(basepath);
rippfile = fullfile(basepath, [basename, '.ripp.mat']);

% load
if isempty(ripp)
    if exist(rippfile, 'file')
        load(rippfile, 'ripp')
    else
        error('ripp file missing')
    end
end

% load session vars
varsFile = ["session"; "spikes.cellinfo"; "spktimes"];
varsName = ["session"; "spikes"; "spktimes"];
v = getSessionVars('basepaths', {basepath}, 'varsFile', varsFile,...
    'varsName', varsName);

% params
fsSpk       = v.session.extracellular.sr;
fs          = ripp.info.fs;
epochs      = ripp.epochs;
nepochs     = size(epochs, 1);
nbinsMap    = size(ripp.maps.freq, 2);
durWin      = [-75 75] / 1000;
recWin      = ripp.info.recWin;
if recWin(2) > floor(v.session.general.duration)
    recWin(2) = floor(v.session.general.duration);
end

% timebins 
tbins = [0, 1 / fs : 1 / fs : diff(recWin)];
tbins = tbins(:) + recWin(1);

% initialize
ripp.spks = [];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create non-ripple epochs for baseline comparison
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get recording epochs without ripples, including confidance margin
nonRipp = SubtractIntervals([recWin(1), recWin(2)], epochs + durWin);        

% make sure these epochs are long enough
durIdx = diff(nonRipp') < sum(abs(durWin));
nonRipp(durIdx, :) = [];

% create timestamps without ripples
randTimes = Restrict([tbins], nonRipp);

% make sure the number of non-ripple epochs is at least as nepochs. If not
% that most definitely means there was a problem with the detection
if length(randTimes) < nepochs
    error('Check ripple detection')
end

% select a random subset of non-ripples epochs
randIdx = randperm(length(randTimes), nepochs);

% take the center of non-ripple epochs as the "peak position"
randPos = sort(randTimes(randIdx));

% make sure this worked
if ~all(Restrict(randPos, nonRipp)) || any(InIntervals(nonRipp, epochs))
    error('finding non-ripples epochs failed')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic 
if ~isempty(v.spktimes)      
    fprintf('Getting MU spikes in ripples... ')

    muSpks = sort(vertcat(v.spktimes{:})) / fsSpk;
    if fullAnalysisFlag
        % find all spikes from all tetrodes that occured during each
        % ripple. save their absolute timestamp (relative to the recording)
        % and relative timestamps (to the ripples start)
        for iepoch = 1 : nepochs
            ripp.spks.mu.rippAbs{iepoch} = muSpks(muSpks < epochs(iepoch, 2) &...
                muSpks > epochs(iepoch, 1));
            ripp.spks.mu.rippRel{iepoch} = ripp.spks.mu.rippAbs{iepoch} - epochs(iepoch, 1);
        end

        % find the mu spike rate per ripple by dividing the total number of
        % spikes with the ripple duration
        nspksRipp = cellfun(@length, ripp.spks.mu.rippAbs, 'uni', true);
        ripp.spks.mu.rate = nspksRipp ./ ripp.dur';
    end

    % count spikes in bins corresponding to the lfp sampling frequency
    spkRate = histcounts(muSpks, [tbins']);

    % create map of spiking activity during ripples
    [r, i] = Sync([tbins(2 : end) spkRate'], ripp.peakPos,...
        'durations', ripp.maps.durWin);
    ripp.spks.mu.rippMap = SyncMap(r, i, 'durations', ripp.maps.durWin,...
        'nbins', nbinsMap, 'smooth', 0);

    % create map of spiking activity during random times
    [r, i] = Sync([tbins(2 : end) spkRate'], randPos,...
        'durations', ripp.maps.durWin);
    ripp.spks.mu.randMap = SyncMap(r, i, 'durations', ripp.maps.durWin,...
        'nbins', nbinsMap, 'smooth', 0);

    fprintf('done in %.2f\n', toc)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% su
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isempty(v.spikes)  
    tic
    fprintf('Getting SU spikes in ripples... ')
    nunits = length(v.spikes.times);

    if fullAnalysisFlag
        spks.su = bz_getRipSpikes('basepath', basepath,...
            'events', epochs, 'spikes', v.spikes, 'saveMat', false);

        ripp.spks.su.rippMap = zeros(nunits, nepochs, nbinsMap);
        ripp.spks.su.randMap = zeros(nunits, nepochs, nbinsMap);
        for iunit = 1 : nunits
            nspksRipp = cellfun(@length, spks.su.UnitEventAbs(iunit, :),...
                'uni', true);
            spks.su.rate(iunit, :) = nspksRipp ./ ripp.dur';
        end
    end

    for iunit = 1 : nunits

        % count spikes in sampling frequency of lfp
        spkRate = histcounts(v.spikes.times{iunit}, [tbins']);

        % create map of firing rate during ripples
        [r, i] = Sync([tbins(2 : end) spkRate'], ripp.peakPos,...
            'durations', ripp.maps.durWin);
        ripp.spks.su.rippMap(iunit, :, :) = SyncMap(r, i, 'durations',...
            ripp.maps.durWin, 'nbins', nbinsMap, 'smooth', 0);

        % create map of firing rate during random times
        [r, i] = Sync([tbins(2 : end) spkRate'], randPos,...
            'durations', ripp.maps.durWin);
        ripp.spks.su.randMap(iunit, :, :) = SyncMap(r, i, 'durations',...
            ripp.maps.durWin, 'nbins', nbinsMap, 'smooth', 0);
    end

    fprintf('done in %.2f\n', toc)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save and graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if saveVar
    save(rippfile, 'ripp', '-v7.3')
end

if graphics 
    plot_rippleSpks(ripp, 'basepath', basepath, 'saveFig', true)
end

end

% EOF