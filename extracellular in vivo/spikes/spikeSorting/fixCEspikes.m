function spikes = getSpikes(varargin)
% bz_getSpikes - Get spike timestamps.
%
% USAGE
%
%    spikes = bz_getSpikes(varargin)
% 
% INPUTS
%
%    spikeGroups     -vector subset of shank IDs to load (Default: all)
%    region          -string region ID to load neurons from specific region
%                     (requires sessionInfo file or units->structures in xml)
%    UID             -vector subset of UID's to load 
%    basepath        -path to recording (where .dat/.clu/etc files are)
%    getWaveforms    -logical (default=true) to load mean of raw waveform data
%    forceL          -logical (default=false) to force loading from
%                     res/clu/spk files
%    saveMat         -logical (default=false) to save in buzcode format
%    noPrompts       -logical (default=false) to supress any user prompts
%    
% OUTPUTS
%
%    spikes - cellinfo struct with the following fields
%          .sessionName    -name of recording file
%          .UID            -unique identifier for each neuron in a recording
%          .times          -cell array of timestamps (seconds) for each neuron
%          .spindices      -sorted vector of [spiketime UID], useful for 
%                           input to some functions and plotting rasters
%          .region         -region ID for each neuron (especially important large scale, high density probes)
%          .shankID        -shank ID that each neuron was recorded on
%          .maxWaveformCh  -channel # with largest amplitude spike for each neuron
%          .rawWaveform    -average waveform on maxWaveformCh (from raw .dat)
%          .cluID          -cluster ID, NOT UNIQUE ACROSS SHANKS
%           
% NOTES
%
% This function can be used in several ways to load spiking data.
% Specifically, it loads spiketimes for individual neurons and other
% sessionInfodata that describes each neuron.  Spiketimes can be loaded using the
% UID(1-N), the shank the neuron was on, or the region it was recorded in.
% The default behavior is to load all spikes in a recording. The .shankID
% and .cluID fields can be used to reconstruct the 'units' variable often
% used in FMAToolbox.
% units = [spikes.shankID spikes.cluID];
% 
% 
% first usage recommendation:
% 
%   spikes = bz_getSpikes('saveMat',true); Loads and saves all spiking data
%                                          into buzcode format .cellinfo. struct
% other examples:
%
%   spikes = bz_getSpikes('spikeGroups',1:5); first five shanks
%
%   spikes = bz_getSpikes('region','CA1'); cells tagged as recorded in CA1
%
%   spikes = bz_getSpikes('UID',[1:20]); first twenty neurons
%
%
% written by David Tingley, 2017
% 
% 23 nov 18 LH - added avg and std waveform fields
% 07 feb 20 LH - fixed orientation given one electrode

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addParameter(p, 'basepath', pwd, @ischar);
addParameter(p, 'spikes', [], @isstruct);
addParameter(p, 'force', false, @islogical);
addParameter(p, 'saveVar', true, @islogical);

parse(p, varargin{:})
basepath = p.Results.basepath;
force = p.Results.force;
saveVar = p.Results.saveVar;

cd(basepath)
session = CE_sessionTemplate(pwd, 'viaGUI', false,...
    'force', false, 'saveVar', false);
fs = session.extracellular.sr;
spkgrp = session.extracellular.spikeGroups.channels;
nsamps = ceil(1.6 * 10^-3 * fs);        % 1.6 ms in samples

% load spikes struct
[~, filename] = fileparts(basepath);
spikesname = [filename '.spikes.cellinfo.mat'];
if exist(spikesname, 'file') 
    fprintf('\nloading %s\n', spikesname)
    load(spikesname)
else
    fprintf('\nextracting spike waveforms\n')
end

% check that fields exist
if exist('spikes', 'var')
    if all(isfield(spikes, {'avgwv', 'stdwv'})) && force == false
        fprintf('all fields exist, skipping...')
        return
    end
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% handle files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% find res / clu / fet / spk files here
clufiles = dir([basepath filesep '*.clu*']);  
resfiles = dir([basepath filesep '*.res*']);
spkfiles = dir([basepath filesep '*.spk*']);

% remove *temp*, *autosave*, and *.clu.str files/directories
tempFiles = zeros(length(clufiles), 1);
for i = 1 : length(clufiles) 
    dummy = strsplit(clufiles(i).name, '.'); 
    if ~isempty(findstr('temp', clufiles(i).name)) |...
            ~isempty(findstr('autosave', clufiles(i).name)) |...
            isempty(str2num(dummy{length(dummy)})) |...
        find(contains(dummy, 'clu')) ~= length(dummy)-1  
        tempFiles(i) = 1;
    end
end
clufiles(tempFiles == 1) = [];

tempFiles = zeros(length(resfiles), 1);
for i = 1:length(resfiles)
    if ~isempty(findstr('temp', resfiles(i).name)) |...
            ~isempty(findstr('autosave', resfiles(i).name))
        tempFiles(i) = 1;
    end
end
resfiles(tempFiles == 1) = [];
tempFiles = zeros(length(spkfiles), 1);
for i = 1:length(spkfiles)
    if ~isempty(findstr('temp', spkfiles(i).name)) |...
            ~isempty(findstr('autosave', spkfiles(i).name))
        tempFiles(i) = 1;
    end
end
spkfiles(tempFiles == 1) = [];

if isempty(clufiles) || isempty(resfiles) || isempty(spkfiles)
    error('no files found...')
end

if length(resfiles) ~= length(clufiles) ||...
        length(clufiles) ~= length(spkfiles)
    warning('different number of res / clu / spk / fet files')
end

% ensures load in sequential order
for i = 1 : length(clufiles)
    temp = strsplit(clufiles(i).name, '.');
    shanks(i) = str2num(temp{length(temp)});
end
[shanks, ind] = sort(shanks);
clufiles = clufiles(ind); 
resfiles = resfiles(ind);
spkfiles = spkfiles(ind);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k = 1;
for i = 1 : length(clufiles)     

    temp = strsplit(clufiles(i).name,'.');
    shankID = str2num(temp{length(temp)});
    electrodes = spkgrp{shankID};       % note this is 1-based
    nchans = length(electrodes);
    fprintf('working on shank %d...\n', shankID)
    
    % clu
    cluname = fullfile(clufiles(i).folder, clufiles(i).name);
    fid = fopen(cluname, 'r');
    nclu = fscanf(fid, '%d', 1);
    clu = fscanf(fid, '%f', Inf);
    fclose(fid);
    cluGrp = unique(clu);
    if length(cluGrp) ~= nclu
        warning('cluFile has the wrong number of clusters')
        nclu = length(cluGrp);
    end
    nspks = length(clu);
    
    % res
    resname = fullfile(resfiles(i).folder, resfiles(i).name);
    fid = fopen(resname, 'r');
    res = fscanf(fid, '%d', Inf);
    fclose(fid);
    
    % spk
    spkname = fullfile(spkfiles(i).folder, spkfiles(i).name);
    fid = fopen(spkname, 'r');
    spk = fread(fid, Inf, 'int16');
    spk = reshape(spk, nchans, nsamps, nspks);
    fclose(fid);
    
    % remove noise clusters
    cluGrp(cluGrp == 0) = [];
    cluGrp(cluGrp == 1) = [];  
      
    % arrange data
    for ii = 1 : length(cluGrp)
        cluidx = find(clu == cluGrp(ii));
        spikes.ts{k} = res(cluidx);
        spikes.times{k} = spikes.ts{k} ./ fs;
        spikes.UID(k) = k;
        spikes.cluID(k) = cluGrp(ii);
        spikes.shankID(k) = shankID;
                    
        avgwv = mean(spk(:, :, cluidx), 3);
        stdwv = std(spk(:, :, cluidx), [], 3);
        [~, ampMaxCh] = max(range(avgwv, 2));

        spikes.rawWaveform{k} = avgwv(ampMaxCh, :);
        spikes.rawWaveform_all{k} = avgwv;
        spikes.rawWaveform_std{k} = stdwv(ampMaxCh, :);
        spikes.rawWaveform_std_all{k} = stdwv;
        spikes.filtWaveform{k} = avgwv(ampMaxCh, :);
        spikes.filtWaveform_all{k} = avgwv;
        spikes.filtWaveform_std{k} = stdwv(ampMaxCh, :);
        
        spikes.maxWaveformCh1(k) = spkgrp{i}(ampMaxCh);
        spikes.maxWaveformCh(k) = spkgrp{i}(ampMaxCh) - 1;
            
        k = k + 1;
    end   
end

% save to buzcode format (before exclusions)
if saveVar
    save(spikesname, 'spikes')
end

end

% EOF

