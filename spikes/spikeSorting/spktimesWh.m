function [spktimes, spkch] = spktimesWh(varargin)

% detects spikes from whitened data. based on ks functions, specifically
% isloated_peaks_new.

% INPUT:
%   basepath    string. path to recording folder where temp_wh.dat exists {pwd} 
%   fs          numeric. sampling frequency [hz]{20000}
%   nchans      numeric. number of channels in dat file.
%   spkgrp      array where each cell is the electrodes for a spike group. 
%   chunksize   numeric. size of chunk to load [samples]. note data 
%               is loaded and processed in batces of spkgrps
%   winCalc     numeric. time range to whiten [s]{[0 Inf]}
%   winWh       numeric. time range to use for the cov mat [s]{[0 Inf]}
%   saveWh      logical. save temp_wh.dat {false}
%   saveVar     logical. save output {true}
%   graphics    logical. plot graphics {false}
%   force       logical. perform detection even if spktimes exists {false}
%   createWh    logical. create temp_wh.dat even file exists {true}.
%               currently not working
%
% DEPENDENCIES
%   get_whitening_matrix (ks)
%   my_min (ks)
%   getOr (ks)
%   class2bytes
%
% TO DO LIST:
% 
% 01 nov 20 LH      
% 22 aug 21      allows for empty tetrodes. this is so that rogue
%                electrodes do not contaminate the whitening matrix 
% 15 may 22      getWhiteMat from get_whitening_matrix to add window

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addOptional(p, 'basepath', pwd);
addOptional(p, 'fs', 20000, @isnumeric);
addOptional(p, 'nchans', 35, @isnumeric);
addOptional(p, 'spkgrp', {}, @iscell);
addOptional(p, 'chunksize', 1024 ^ 2 + 64, @isnumeric);
addOptional(p, 'winCalc', [0 Inf], @isnumeric);
addOptional(p, 'winWh', [0 Inf], @isnumeric);
addOptional(p, 'saveWh', true, @islogical);
addOptional(p, 'saveVar', true, @islogical);
addOptional(p, 'graphics', false, @islogical);
addOptional(p, 'force', false, @islogical);
addOptional(p, 'createWh', true, @islogical);

parse(p, varargin{:})
basepath    = p.Results.basepath;
fs          = p.Results.fs;
nchans      = p.Results.nchans;
spkgrp      = p.Results.spkgrp;
chunksize   = p.Results.chunksize;
winCalc     = p.Results.winCalc;
winWh       = p.Results.winWh;
saveWh      = p.Results.saveWh;
saveVar     = p.Results.saveVar;
graphics    = p.Results.graphics;
force       = p.Results.force;
createWh    = p.Results.createWh;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% adjust spkgrps for empty tetrodes
[~, chIdx] = sort([spkgrp{:}]);
k = 1;
for igrp = 1 : length(spkgrp)
    for ich = 1 : length(spkgrp{igrp})
        spkgrp_ch{igrp}(ich) = chIdx(k);
        k = k + 1;
    end
end

% constants
loc_range = [6 max(size(cell2nanmat(spkgrp, 2), 1))];  % for running minimum (dim 1 - sample; dim 2 - channel] 
scaleproc = 200;    % conversion back from unit variance to voltage
thr = -8;           % spike threshold in standard deviations {-6}
dt = 8;             % dead time between spikes [samples]

% whitening params from ks
ops = opsKS('basepath', basepath, 'fs', fs, 'nchans', nchans,...
    'spkgrp', spkgrp, 'trange', winCalc);
NT              = chunksize;        % override
ops.NT          = NT;
NchanTOT        = ops.NchanTOT; 
bytes           = dir(ops.fbinary);
bytes           = bytes.bytes;
nTimepoints     = floor(bytes / NchanTOT / 2); 
ops.tstart      = ceil(ops.trange(1) * ops.fs); 
ops.tend        = min(nTimepoints, ceil(ops.trange(2) * ops.fs)); 
ops.sampsToRead = ops.tend-ops.tstart; 
ops.twind       = ops.tstart * NchanTOT*2; 
Nbatch          = ceil(ops.sampsToRead / NT); 
ops.Nbatch      = Nbatch;
load(ops.chanMap)
ops.igood       = true(size(chanMap));
ops.Nchan       = numel(kcoords); 
ops.Nfilt       = getOr(ops, 'nfilt_factor', 4) * ops.Nchan; 
rez.ops         = ops; 
rez.xc          = xcoords; % for historical reasons, make all these copies of the channel coordinates
rez.yc          = ycoords;
rez.xcoords     = xcoords;
rez.ycoords     = ycoords;
rez.ops.kcoords = kcoords;
NTbuff          = NT + 3 * ops.ntbuff; % we need buffers on both sides for filtering
rez.ops.NTbuff  = NTbuff;
rez.ops.chanMap = chanMap;

% check if spktimes already exists
[~, basename] = fileparts(basepath);
varname = fullfile(basepath, [basename '.spktimes.mat']);
if exist(varname) && ~force
    load(varname)
    return
end

% check if temp_wh already exists
if ~exist(ops.fproc, 'file')
    createWh = true;
end

fprintf('/n/ndetecting spikes in %s/n', basename)
tic;

if createWh

    % weights to combine batches at the edge
    w_edge          = linspace(0, 1, ops.ntbuff)';
    ntb             = ops.ntbuff;
    datr_prev       = gpuArray.zeros(ntb, ops.Nchan, 'single');

    % get a rotation matrix (Nchan by Nchan) which whitens the zero-timelag
    % covariance of the data
    if winWh(2) > nTimepoints / fs
        winWh(2) = floor(nTimepoints / fs);
    end
    Wrot = getWhiteMat(rez, winWh);

    fid = fopen(ops.fbinary, 'r'); % open for reading raw data
    if saveWh
        fidW = fopen(ops.fproc, 'w+'); % open for writing processed data
    end

else
    fidW = fopen(ops.fproc, 'r'); % open for writing processed data

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% detection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize
spktimes = cell(1, length(spkgrp));
spkch = cell(1, length(spkgrp));
spkamp = cell(1, length(spkgrp));
for ibatch = 1 : Nbatch          
     
    % print progress
    if ibatch ~= 1
        fprintf(repmat('\b', 1, length(txt)))
    end
    txt = sprintf('working on chunk %d / %d', ibatch, Nbatch);
    fprintf(txt)
    
    % create or load whitened data
    if createWh
        offset = max(0, ops.twind + 2 * NchanTOT * (NT * (ibatch - 1) - ntb));
        fseek(fid, offset, 'bof');
        buff = fread(fid, [NchanTOT NTbuff], '*int16');
        nsampcurr = size(buff,2);
        if nsampcurr<NTbuff
            buff(:, nsampcurr + 1 : NTbuff) = repmat(buff(:, nsampcurr), 1, NTbuff - nsampcurr);
        end
        if offset == 0
            bpad = repmat(buff(:, 1), 1, ntb);
            buff = cat(2, bpad, buff(:, 1 : NTbuff-ntb));
        end
        wh = gpufilter(buff, ops, chanMap);
        wh(ntb + [1 : ntb], :) = w_edge .* wh(ntb + [1 : ntb], :) +...
            (1 - w_edge) .* datr_prev;
        datr_prev = wh(ntb + NT + [1:ops.ntbuff], :);
        wh    = wh(ntb + (1 : NT), :);
        wh    = wh * Wrot; % whiten the data and scale by 200 for int16 range
        if saveWh
            datcpu  = gather(int16(wh')); % convert to int16, and gather on the CPU side
            fwrite(fidW, datcpu, 'int16'); % write this batch to binary file
        end

    else
        offset = max(0, ops.twind + 2 * NchanTOT * (NT * (ibatch - 1)));
        fseek(fidW, offset, 'bof');
        wh = fread(fidW, [NchanTOT (NTbuff - ntb * 3)], '*int16')';
        wh = wh(:, [spkgrp{:}]);
        wh    = wh(ntb + (1 : NT), [spkgrp{:}]);
        wh = wh / scaleproc;

    end

    wh = wh / scaleproc;   
    
    % find spikes on each tetrode 
    for igrp = 1 : length(spkgrp_ch)
        
        % skip empty tetrodes
        if isempty(spkgrp_ch{igrp})
            continue
        end
        
        % moving minimum across channels and samples 
        smin = my_min(wh(:, spkgrp_ch{igrp}), loc_range, [1, 2]);
        
        % peaks are samples that achieve this local minimum AND have
        % negativities less than a preset threshold
        crossings = wh(:, spkgrp_ch{igrp}) < smin + 1e-5 &...
            wh(:, spkgrp_ch{igrp}) < thr;
        
        % find the non-zero peaks, and take their amplitudes.
        [samp, ch, amp] = find(crossings .* wh(:, spkgrp_ch{igrp})); 
        [samp, ia] = sort(gather(samp + NT * (ibatch - 1)));
        ch = gather(ch(ia));
        amp = gather(amp(ia));
                
        % if two spikes are closer than dt samples, keep larger one. this
        % complemants the local minimum in cases where the raw peak is wide
        % such that more than one sample reaches the local minimum.
        while any(diff(samp) < dt)
            idx = find(diff(samp) <= dt);
            idx2 = amp(idx) > amp(idx + 1);
            idx3 = sort([idx(idx2); idx(~idx2) + 1]);
            amp(idx3) = [];
            samp(idx3) = [];
            ch(idx3) = [];
        end        
        spktimes{igrp} = [spktimes{igrp}; samp]; 
        spkch{igrp} = [spkch{igrp}; spkgrp_ch{igrp}(ch)'];
        if graphics
            spkamp{igrp} = [spkamp{igrp}; amp];
        end
    end
end

% save variable
if saveVar
    save(varname, 'spktimes', 'spkch')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if graphics
   
    % ---------------------------------------------------------------------
    % inspect detection 
       
    % raw
    buff = single(buff');
    buff = buff(ntb + (1:NT),:);   
     
    % params
    xLimit = [1 : length(buff)]';
    clr = ['kbgroym'];
    grp = 3;            % selected grp        
    
    figure
    ax1 = subplot(3, 1, 1);
    yOffset = gather(median(range(buff)));
    hold on
    for igrp = length(spkgrp{grp}) : -1 : 1
        plot(xLimit / fs, buff(:, spkgrp{grp}(igrp))...
            + yOffset * (igrp - 1), clr(igrp))
    end
    set(gca, 'TickLength', [0 0])
    title('Raw')
    
    % wh
    ax2 = subplot(3, 1, 2);
    yOffset = gather(max(range(wh)));
    hold on
    for igrp = length(spkgrp{grp}) : -1 : 1
        plot(xLimit / fs, wh(:, spkgrp{grp}(igrp))...
            + yOffset * (igrp - 1), clr(igrp))
        scatter((spktimes{grp}(spkch{grp} == spkgrp{grp}(igrp)) - NT * (igrp - 1)) / fs,...
            spkamp{grp}(spkch{grp} == spkgrp{grp}(igrp))...
            + yOffset * (igrp - 1), '*')
        yline(thr + yOffset * (igrp - 1), '--r');
    end
    set(gca, 'TickLength', [0 0])
    title('Whitened')
    
    % min
    smin = my_min(wh(:, spkgrp{grp}), loc_range, [1, 2]);
    ax3 = subplot(3, 1, 3);
    yOffset = gather(max(range(smin)));
    hold on
    for igrp = length(spkgrp{grp}) : -1 : 1
        plot(xLimit / fs, smin(:, igrp)...
            + yOffset * (igrp - 1), clr(igrp))
    end
    set(gca, 'TickLength', [0 0])
    xlabel('Time [s]')
    title('Local Minimum')
    linkaxes([ax1, ax2, ax3], 'x')
    
    % ---------------------------------------------------------------------
    % isi histogram
    figure
    binEdges = [0 : dt : 1000];      % 0 - 50 ms
    histogram(diff(samp), binEdges, 'EdgeColor', 'none', 'FaceColor', 'k')
    xlim([binEdges(1) binEdges(end)])
    xlabel('ISI [samples]')
    ylabel('Counts')    
end

% close files
fclose(fid);
if saveWh
    fclose(fidW);
end

fprintf('\nthat took %.2f minutes\n', toc / 60)

end

% EOF