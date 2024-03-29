function LFPfromDat(varargin)

% performs lowpass filter on wideband data (basename.dat file), subsamples
% the filtered data and saves as a new flat binary. filtering is performed
% by FFT-based convolution with the sinc kernel (IOSR.DSP.SINCFILTER).
% based on bz_LFPdromDat with the following differences: (1) no
% dependencies on sessionInfo, (2) bz handles the remainder separately s.t.
% there is a continuity problem at the last chunk, (3) feeds sincFilter all
% channels at once (not much difference because sincFilter loops through
% channels anyway), (4) slightly faster (180 vs 230 s for 90 m recording),
% (5) no annoying waitbar.
%  
% because the way downsampling occurs, the final duration may not be
% accurate if the input fs is not a round number. one way arround
% this is to define the output fs s.t. the ratio output/input is round.
% however, this is not possible for many inputs (including 24414.06). do
% not see a way around this inaccuracy (~50 ms for 90 m recording). 
% 
% despite the description in IOSR.DSP.SINCFILTER, the output is slightly
% different if the cutoff frequency is [0 450] instead of [450]. thus only
% a scalar (i.e. cutoff for low-pass) is allowed here. 
% 
% INPUT
%   basepath    string. path to load filename and save output {pwd}
%   precision   char. sample precision of dat file {'int16'} 
%   clip        mat n x 2 indicating samples to diregard from chunks.
%               for example: clip = [0 50; 700 Inf] will remove the first
%               50 samples and all samples between 700 and n
%   fsIn        numeric [Hz]. sampling frequency of dat file {20000}
%   fsOut       numeric [Hz]. sampling frequency of lfp file {1250}
%   nchans      numeric. number of channels in dat file {16}
%   cf          scalar. cutoff frequency {[450]}.
%   chunksize   size of data to load at once [samples]{5e6}. 
%               if empty will load entire file (be careful!).
%               for 35 channels in int16, 5e6 samples = 350 MB.
%   force       logical. force load even if file exists {false}
%   
% DEPENDENCIES
%   IOSR.DSP.SINCFILTER
%   class2bytes
% 
% 11 aug 20 LH
%
% TO DO LIST
%   # embed sinc filter (to avoid repetitive calculation of kernel)       


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic

p = inputParser;
addOptional(p, 'basepath', pwd);
addOptional(p, 'precision', 'int16', @ischar);
addOptional(p, 'nchans', 35, @isnumeric);
addOptional(p, 'clip', [], @isnumeric);
addOptional(p, 'fsIn', 20000, @isnumeric);
addOptional(p, 'fsOut', 1250, @isnumeric);
addOptional(p, 'cf', 450, @isnumeric);
addOptional(p, 'chunksize', 5e6, @isnumeric);
addOptional(p, 'force', false, @islogical);

parse(p, varargin{:})
basepath = p.Results.basepath;
precision = p.Results.precision;
nchans = p.Results.nchans;
clip = p.Results.clip;
fsIn = p.Results.fsIn;
fsOut = p.Results.fsOut;
cf = p.Results.cf;
chunksize = p.Results.chunksize;
force = p.Results.force;

% suppress integer warning (for tdt fs)
warning('off', 'MATLAB:colon:nonIntegerIndex') 

% import sync filter toolbox
import iosr.dsp.*

% size of one data point in bytes
nbytes = class2bytes(precision);

filtRatio = cf / (fsIn / 2);
%  from SINCFILTER: If cf is a scalar, then cf specifies the low-pass
%  cutoff frequency. If cf is a two-element vector, then cf specifies the
%  band-pass interval. cf must be 0.0 < cf < 1.0, with 1.0 corresponding to
%  half the sample rate. 

fsRatio = (fsIn / fsOut);
if cf > fsOut / 2
    warning('low pass cutoff beyond nyquist')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% handle files and chunks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% handle files
[~, basename] = fileparts(basepath);
fdat = fullfile(basepath,[basename,'.dat']);
flfp = fullfile(basepath,[basename,'.lfp']);

% check that basename.dat exists
if ~exist(fdat, 'file')
    error('%s does not exist', fdat)
end
datinfo = dir(fdat);

% check if basename.lfp exists
if exist(flfp, 'file') && ~force
  fprintf('%s exists, returning...\n\n', flfp)
  return
end

% Set chunk and buffer size as even multiple of fsRatio
if mod(chunksize, fsRatio) ~= 0
    chunksize = round(chunksize + fsRatio - mod(chunksize, fsRatio));
end

ntbuff = 525;  % default filter size in iosr toolbox
if mod(ntbuff, fsRatio) ~= 0
    ntbuff = round(ntbuff + fsRatio - mod(ntbuff, fsRatio));
end

% partition into chunks
nsamps = datinfo.bytes / nbytes / nchans;
chunks = n2chunks('n', nsamps, 'chunksize', chunksize, 'clip', clip,...
    'overlap', ntbuff);
nchunks = size(chunks, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\ncreating %s\n', flfp)

% memory map to original file
m = memmapfile(fdat, 'Format', {precision [nchans nsamps] 'mapped'});
raw = m.data;

fid = fopen(fdat, 'r');
fidOut = fopen(flfp, 'a');
for i = 1 : nchunks
    
    % print progress
    if i ~= 1
        fprintf(repmat('\b', 1, length(txt)))
    end
    txt = sprintf('working on chunk %d / %d', i, nchunks);
    fprintf(txt)
    
    % load chunk
    d = raw.mapped(:, chunks(i, 1) : chunks(i, 2));
    d = double(d);
    
    % filter
    filtered = [iosr.dsp.sincFilter(d', filtRatio)]';
    
    % downsample
    if i == 1
        dd = int16(real(filtered(:, fsRatio : fsRatio :...
            length(filtered) - ntbuff)));
    else
        dd = int16(real(filtered(:, ntbuff + fsRatio : fsRatio :...
            length(filtered) - ntbuff)));
    end

    fwrite(fidOut, dd(:), 'int16'); 
end

clear raw
clear m
fclose(fid);
fclose(fidOut);

fprintf('\nLFP file created. that took %.2f minutes\n', toc / 60)

% restore integer warning (for tdt fs)
warning('on', 'MATLAB:colon:nonIntegerIndex') 

end

% EOF




% EOF