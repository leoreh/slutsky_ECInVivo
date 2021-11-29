function burstiness = burstMetrics(spikes, varargin)

%  calculate metrics of burstiness from ACG. based on calc_ACG_metrics from
%  cell explorer. metrics include: busrtIndex_Royer2012, Mizuseki2012,
%  Doublets. calculates two ACGs: narrow (100ms, 0.5ms bins) and wide (1s,
%  1ms bins)


% INPUT:
%   spikes          struct (see getSpikes)
%   sunits          numeric vec. indices of selected units for calculation
%                   {[]}.
%   basepath        path to recording
%   graphics        logical. plot graphics {true} or not (false)
%   saveFig         logical. save figure {1} or not (0)
%   saveVar         logical. save variables (update spikes and save su)
%   session         struct. session info (see CE_sessionTemplate)
%
% OUTPUT:           

%
% DEPENDENCIES:
%   CCG             
%
% TO DO LIST:
%
% 24 nov 21 LH      

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addOptional(p, 'spikes', []);
addOptional(p, 'sunits', []);
addOptional(p, 'basepath', pwd, @ischar);
addOptional(p, 'graphics', true, @islogical);
addOptional(p, 'saveFig', true, @islogical);
addOptional(p, 'saveVar', true, @islogical);

parse(p, varargin{:})
spikes      = p.Results.spikes;
basepath    = p.Results.basepath;
graphics    = p.Results.graphics;
saveFig     = p.Results.saveFig;
saveVar     = p.Results.saveVar;

% load spikes if empty
if isempty(spikes)
    [~, filename] = fileparts(basepath);
    spkname = [filename '.spikes.cellinfo.mat'];
    if exist(spkname, 'file')
        load(spkname)
    else
        error('%s not found', spkname)
    end
end

% make sure spikes has required fields
if ~all(isfield(spikes, {'shankID', 'cluID', 'times'}))
    error('spikes missing required fields')
end

% selected untis
if isemtpy(sunits)
    sunits = 1 : length(spikes.times);
end

% load session info
[~, basename] = fileparts(basepath);
sessionName = [basename, '.session.mat'];
if ~exist(sessionName, 'file')
    session = CE_sessionTemplate(pwd, 'viaGUI', false,...
        'force', true, 'saveVar', true);
else
    load(sessionName)
end

% params
fs = session.extracellular.sr;
spkgrp = session.extracellular.spikeGroups.channels;
ngrp = length(unique(spikes.shankID));  % only tetrodes with units
nunits = length(spikes.times);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calc ACGs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize
bins_wide = 500;
binsize_wide = 0.001;
duration_wide = 1; 
acg_wide = zeros(duration_wide / binsize_wide + 1, length(sunits));

bins_narrow = 100;
binsize_narrow = 0.0005;
duration_narrow = 0.1; 
acg_narrow = zeros(duration_narrow / binsize_narrow + 1, length(sunits));

for iunit = sunits
    acg_wide(:, iunit) = CCG(spikes.times{iunit},...
        ones(size(spikes.times{iunit})), 'binSize', binsize_wide,...
        'duration', duration_wide, 'norm', 'rate', 'Fs', 1 / fs);
    
    acg_narrow(:, iunit) = CCG(spikes.times{iunit},...
        ones(size(spikes.times{iunit})), 'binSize', binsize_narrow,...
        'duration', duration_narrow, 'norm', 'rate', 'Fs', 1 / fs);
    
    % doublets: max bin count from 2.5-8ms normalized by the average number
    % of spikes in the 8-11.5ms bins
     burstiness.doublets(iunit) = max(acg_narrow(bins_narrow + 1 + 5 : bins_narrow + 1 + 16, iunit)) /...
         mean(acg_narrow(bins_narrow + 1 + 16 : bins_narrow + 1 + 23, iunit));
    
    % royer 2012: average number of spikes in the 3-5 ms bins divided by the
    % average number of spikes in the 200-300 ms bins
     burstiness.royer2012(iunit) = mean(acg_wide(bins_wide + 1 + 3 : bins_wide + 1 + 5, iunit)) /...
         mean(acg_wide(bins_wide + 1 + 200 : bins_wide + 1 + 300, iunit));
     
    % Lv (firing irregularity): adapted from Shinomoto 2003 and Kobayashi 2019.
     
    % Mizuseki 2011: fraction of spikes with a ISI for following or preceding
    % spikes < 0.006
    bursty = [];
    for jj = 2 : length(spikes{spkExclu}.times{j}) - 1
        bursty(jj) =  any(diff(spikes{spkExclu}.times{j}(jj-1 : jj + 1)) < 0.006);
    end
    cell_metrics.burstIndex_Mizuseki2012(j) = length(find(bursty > 0))/length(bursty); % 
     
end

figure
subplot(3,1,2)
histogram(burstiness.royer2012, 40),xlabel('BurstIndex Royer2012'), ylabel('Count')
subplot(3,1,3)
histogram(burstiness.doublets, 40),xlabel('BurstIndex Doublets'), ylabel('Count')




end

% EOF