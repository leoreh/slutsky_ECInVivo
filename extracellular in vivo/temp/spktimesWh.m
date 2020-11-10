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
%   saveWh      logical. save temp_wh.mat {false}
%   saveVar     logical. save output {true}
%   graphics    logical {false}
%
% DEPENDENCIES
%   get_whitening_matrix (ks)
%   my_min (ks)
%   loadChanMap (ks)
%   get_file_size (ks)
%   class2bytes
%
% TO DO LIST:
% 
% 01 nov 20 LH      

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;

p = inputParser;
addOptional(p, 'basepath', pwd);
addOptional(p, 'fs', 20000, @isnumeric);
addOptional(p, 'nchans', 35, @isnumeric);
addOptional(p, 'spkgrp', {}, @iscell);
addOptional(p, 'chunksize', 2048 ^ 2 + 64, @isnumeric);
addOptional(p, 'saveWh', true, @islogical);
addOptional(p, 'saveVar', true, @islogical);
addOptional(p, 'graphics', false, @islogical);

parse(p, varargin{:})
basepath    = p.Results.basepath;
fs          = p.Results.fs;
nchans      = p.Results.nchans;
spkgrp      = p.Results.spkgrp;
chunksize   = p.Results.chunksize;
saveWh      = p.Results.saveWh;
saveVar     = p.Results.saveVar;
graphics    = p.Results.graphics;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nchansWh = length([spkgrp{:}]);

% constants
loc_range = [6 max(size(cell2nanmat(spkgrp), 1))];  % for running minimum (dim 1 - sample; dim 2 - channel] 
scaleproc = 200;    % conversion back from unit variance to voltage
thr = -4;           % spike threshold in standard deviations {-6}
dt = 8;             % dead time between spikes [samples]

% whitening params from ks
ops = opsKS('basepath', basepath, 'fs', fs, 'nchans', nchans,...
    'spkgrp', spkgrp, 'trange', [0 Inf]);
NT              = ops.NT;
NT              = 1024 ^ 2 + 64;        % override
ops.NT          = NT;
NchanTOT        = ops.NchanTOT; 
bytes           = get_file_size(ops.fbinary); 
nTimepoints     = floor(bytes/NchanTOT/2); 
ops.tstart      = ceil(ops.trange(1) * ops.fs); 
ops.tend        = min(nTimepoints, ceil(ops.trange(2) * ops.fs)); 
ops.sampsToRead = ops.tend-ops.tstart; 
ops.twind       = ops.tstart * NchanTOT*2; 
Nbatch          = ceil(ops.sampsToRead / NT); 
ops.Nbatch      = Nbatch;
[chanMap, xc, yc, kcoords, ~] = loadChanMap(ops.chanMap);
ops.igood       = true(size(chanMap));
ops.Nchan       = numel(chanMap); 
ops.Nfilt       = getOr(ops, 'nfilt_factor', 4) * ops.Nchan; 
rez.ops         = ops; 
rez.xc          = xc; % for historical reasons, make all these copies of the channel coordinates
rez.yc          = yc;
rez.xcoords     = xc;
rez.ycoords     = yc;
rez.ops.chanMap = chanMap;
rez.ops.kcoords = kcoords;
NTbuff          = NT + 3*ops.ntbuff; % we need buffers on both sides for filtering
rez.ops.NTbuff  = NTbuff;
rez.ops.chanMap = chanMap;

% weights to combine batches at the edge
w_edge          = linspace(0, 1, ops.ntbuff)';
ntb             = ops.ntbuff;
datr_prev       = gpuArray.zeros(ntb, ops.Nchan, 'single');

% get a rotation matrix (Nchan by Nchan) which whitens the zero-timelag
% covariance of the data
Wrot = get_whitening_matrix(rez); 

fid = fopen(ops.fbinary, 'r'); % open for reading raw data
if saveWh
    fidW = fopen(ops.fproc,   'w+'); % open for writing processed data
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% detection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize
spktimes = cell(1, length(spkgrp));
spkch = cell(1, length(spkgrp));
spkamp = cell(1, length(spkgrp));
for i = 1 : Nbatch          
     
    % print progress
    if i ~= 1
        fprintf(repmat('\b', 1, length(txt)))
    end
    txt = sprintf('working on chunk %d / %d', i, Nbatch);
    fprintf(txt)
    
    % create whitened data
    offset = max(0, ops.twind + 2*NchanTOT*(NT * (i - 1) - ntb)); 
    fseek(fid, offset, 'bof'); 
    buff = fread(fid, [NchanTOT NTbuff], '*int16'); 
    nsampcurr = size(buff,2); 
    if nsampcurr<NTbuff
        buff(:, nsampcurr+1:NTbuff) = repmat(buff(:,nsampcurr), 1, NTbuff-nsampcurr); 
    end
    if offset==0
        bpad = repmat(buff(:,1), 1, ntb);
        buff = cat(2, bpad, buff(:, 1:NTbuff-ntb)); 
    end    
    wh    = gpufilter(buff, ops, chanMap);    
    wh(ntb + [1:ntb], :) = w_edge .* wh(ntb + [1:ntb], :) +...
        (1 - w_edge) .* datr_prev;
    datr_prev = wh(ntb +NT + [1:ops.ntbuff], :);
    wh    = wh(ntb + (1:NT),:);   
    wh    = wh * Wrot; % whiten the data and scale by 200 for int16 range    
    if saveWh
        datcpu  = gather(int16(wh')); % convert to int16, and gather on the CPU side
        fwrite(fidW, datcpu, 'int16'); % write this batch to binary file
    end      
    wh = wh / scaleproc;   
    
    % find spikes on each tetrode 
    for ii = 1 : length(spkgrp)

        % moving minimum across channels and samples 
        smin = my_min(wh(:, spkgrp{ii}), loc_range, [1, 2]);
        
        % peaks are samples that achieve this local minimum AND have
        % negativities less than a preset threshold
        crossings = wh(:, spkgrp{ii}) < smin + 1e-5 &...
            wh(:, spkgrp{ii}) < thr;
        
        % find the non-zero peaks, and take their amplitudes.
        [samp, ch, amp] = find(crossings .* wh(:, spkgrp{ii})); 
        [samp, ia] = sort(gather(samp + NT * (i - 1)));
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
        spktimes{ii} = [spktimes{ii}; samp]; 
        spkch{ii} = [spkch{ii}; spkgrp{ii}(ch)'];
        if graphics
            spkamp{ii} = [spkamp{ii}; amp];
        end
    end
end

% save variable
if saveVar
    [~, basename] = fileparts(basepath);
    save(fullfile(basepath, [basename '.spktimes.mat']),...
        'spktimes', 'spkch')
end

fprintf('\nthat took %.2f minutes\n', toc / 60)

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
    for i = length(spkgrp{grp}) : -1 : 1
        plot(xLimit / fs, buff(:, spkgrp{grp}(i))...
            + yOffset * (i - 1), clr(i))
    end
    set(gca, 'TickLength', [0 0])
    title('Raw')
    
    % wh
    ax2 = subplot(3, 1, 2);
    yOffset = gather(max(range(wh)));
    hold on
    for i = length(spkgrp{grp}) : -1 : 1
        plot(xLimit / fs, wh(:, spkgrp{grp}(i))...
            + yOffset * (i - 1), clr(i))
        scatter((spktimes{grp}(spkch{grp} == spkgrp{grp}(i)) - NT * (i - 1)) / fs,...
            spkamp{grp}(spkch{grp} == spkgrp{grp}(i))...
            + yOffset * (i - 1), '*')
        yline(thr + yOffset * (i - 1), '--r');
    end
    set(gca, 'TickLength', [0 0])
    title('Whitened')
    
    % min
    smin = my_min(wh(:, spkgrp{grp}), loc_range, [1, 2]);
    ax3 = subplot(3, 1, 3);
    yOffset = gather(max(range(smin)));
    hold on
    for i = length(spkgrp{grp}) : -1 : 1
        plot(xLimit / fs, smin(:, i)...
            + yOffset * (i - 1), clr(i))
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

end