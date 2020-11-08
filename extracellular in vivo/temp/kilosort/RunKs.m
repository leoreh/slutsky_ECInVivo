function rez = runKS(varargin)

% wrapper for running kilosort. loads configuration, runs the algorithm and
% arranges the output in phy or neurosuite format.
% 
% INPUT:
%   basepath    string. path to recording folder {pwd}. if multiple dat
%               files exist only the first will be processed
%   procpath    string. path to where processing will occur.
%               should be a fast ssd
%   fs          numeric. sampling frequency [hz]{20000}
%   nchans      numeric. number of channels in dat file.
%   spkgrp      array where each cell is the electrodes for a spike group. 
%   trange      numeric. time range to sort [s]{[0 Inf]}
%   saveFinal   logical. save final rez file {false}
%   viaGui      logical. run ks via gui (true) or script {false}
%   outFormat   string. output results to 'phy' or neurosuite {'ns'} format
%
% DEPENDENCIES
%   npy-matlab
%   kilosort 2
%   opsKS
%   ks2ns
%
% TO DO LIST:
%   # currently, only the gui can be used in cases where one channel of a
%   tetrode is bad (done - 11 jun 20)
% 
% 22 may 20 LH      updates:
% 16 jun 20             outFormat
% 31 oct 20             updated to KS2.5
% 08 nov 20             split opsKS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;

p = inputParser;
addOptional(p, 'basepath', pwd);
addOptional(p, 'procpath', '');
addOptional(p, 'fs', 20000, @isnumeric);
addOptional(p, 'nchans', 35, @isnumeric);
addOptional(p, 'spkgrp', {}, @iscell);
addOptional(p, 'trange', [0 Inf], @isnumeric);
addOptional(p, 'saveFinal', false, @islogical);
addOptional(p, 'viaGui', false, @islogical);
addOptional(p, 'outFormat', 'ns', @ischar);

parse(p, varargin{:})
basepath    = p.Results.basepath;
procpath    = p.Results.procpath;
fs          = p.Results.fs;
nchans      = p.Results.nchans;
spkgrp      = p.Results.spkgrp;
trange      = p.Results.trange;
saveFinal   = p.Results.saveFinal;
viaGui      = p.Results.viaGui;
outFormat   = p.Results.outFormat;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load configuration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ops = opsKS('basepath', basepath, 'fs', fs, 'nchans', nchans,...
    'spkgrp', spkgrp, 'trange', [0 Inf]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run ks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Looking for data inside %s \n', basepath)

if viaGui
    kilosort
else
    % preprocess data to create temp_wh.dat
    rez = preprocessDataSub(ops);
    
    % NEW STEP TO DO DATA REGISTRATION
    rez = datashift2(rez, 1); % last input is for shifting data
    
    % time-reordering as a function of drift
    rez = clusterSingleBatches(rez);
    
    % intermediate saving because the rest can be resumed after loading rez
    save(fullfile(basepath, 'rez.mat'), 'rez', '-v7.3');
    
    % main tracking and template matching algorithm
    iseed = 1; % ORDER OF BATCHES IS NOW RANDOM, controlled by random number generator
    rez = learnAndSolve8b(rez, iseed);
    
    % OPTIONAL: remove double-counted spikes - solves issue in which individual spikes are assigned to multiple templates.
    % See issue 29: https://github.com/MouseLand/Kilosort/issues/29
    rez = remove_ks2_duplicate_spikes(rez);
    
    % final merges
    rez = find_merges(rez, 1);
    
    % final splits by SVD
    rez = splitAllClusters(rez, 1);
    
    % final splits by amplitudes (from v2.0)
    % rez = splitAllClusters(rez, 0);
    
    % decide on cutoff
    rez = set_cutoff(rez);
    % eliminate widely spread waveforms (likely noise)
    rez.good = get_good_units(rez);
       
    fprintf('found %d good units \n', sum(rez.good > 0))
    
    % save final results
    if saveFinal
        % discard features in final rez file (too slow to save)
        rez.cProj = [];
        rez.cProjPC = [];
        rez.st2 = [];

        % final time sorting of spikes, for apps that use st3 directly
        [~, isort]   = sortrows(rez.st3);
        rez.st3      = rez.st3(isort, :);
        
        % Ensure all GPU arrays are transferred to CPU side before saving to .mat
        rez_fields = fieldnames(rez);
        for i = 1:numel(rez_fields)
            field_name = rez_fields{i};
            if(isa(rez.(field_name), 'gpuArray'))
                rez.(field_name) = gather(rez.(field_name));
            end
        end  
        fprintf('\nsaving final rez file\n')
        save(fullfile(basepath, 'rez.mat'), 'rez', '-v7.3');
    end
    
    % write for manual curation
    switch outFormat
        case 'phy'
            fprintf('\nSaving results to Phy \n')
            % note this function delets all .npy in path
            rezToPhy(rez, basepath); 
            % save channel_shanks file
            writeNPY(kcoords(~isnan(kcoords)), fullfile(basepath, 'channel_shanks.npy'));
        case 'ns'
            fprintf('\nSaving results to Neurosuite format \n')
            ks2ns(rez)
    end
    
    % clean temp_wh file
    delete(fullfile(basepath, 'temp_wh.dat'))
end



