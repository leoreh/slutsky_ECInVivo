function acc = EMGfromACC(varargin)

% get 3-axis accelerometer signal from dat file and return the acceleration
% amplitude. based on
% http://intantech.com/files/Intan_RHD2000_accelerometer_calibration.pdf,
% Dhawale et al., eLife, 2017, and bz_EMGFromLFP. Conversion from g to
% volts is done by substracting the zero-g bias [V] and dividing with the
% accelerometer sensitivity [V / g]. these parameters are based on
% measurments done 27.04.20. After conversion to Volts, the 3-axis
% acceleration data is L2 normalized to yeild a magnitude heuristic. The
% broadband spectrogram of the magnitude is calculated and the pc1 (smoothed and
% normalized 0-1) is taken as the output (EMG equivalent).
%
% INPUT:
%   basepath    string. path to .dat and .npy files {pwd}.
%               if multiple dat files exist than concat must be true
%   fname       string. name of dat file. if empty and more than one dat in
%               path, will be extracted from basepath
%   chunksize   size of data to load at once [samples]{5e6}.
%               if empty will load entire file (be careful!).
%               for 35 channels in int16, 5e6 samples = 350 MB.
%   precision   char. sample precision {'int16'}  of dat file
%   clip        mat n x 2 indicating samples to diregard from chunks.
%               for example: clip = [0 50; 700 Inf] will remove the first
%               50 samples and all samples between 700 and n
%   nchans      numeric. number of channels in dat file {35}.
%   ch          3x1 vec. channel index to acceleration in x,
%               y, z (order does not matter so long corresponds to
%               sensitivity and gbias).
%   sensitivity 3x1 vec. sensitivity data to acceleration in x,
%               y. order must correspond to ch.
%   gbias       3x1 vec. zero-g bias data to acceleration in x,
%               y. order must correspond to ch.
%   smf         numeric. 
%   force       logical. reload data even if .mat file exists {false}
%   fsOut       numeric. requested sampling frequency {1250}
%   fsIn        numeric. original sampling frequency {20000}
%   graphics    logical. plot figure {1}
%   saveVar     logical. save variable {1}
%
% OUTPUT
%   acc             struct with the following fields
%       tstamps     time stamps for pc1 [seconds]
%       data        pc1 of spectrogram 
%       fs          sampling frequency
%       fs_orig     original sampling frequency
%
% CALLS:
%   class2bytes
%   specBand
%   binary2epochs
%   sepBimodel
%   bz_BasenameFromBasepath
%   n2chunks
%   bz_NormToRange
%
% TO DO LIST:
%   # compare with emgFromLfp (done)
%   # linear envelop (smoothing and filtering) (not needed)
%   # use pc1 of spectrogram instead of power in 2-5 Hz band (done)
%   # fix sepBimodel (see ClusterStates_GetMetrics for good example)
%   # get data from lfp instead of downsampling dat
%   # update graphics
%   # get signal from LFP file
%
% 09 apr 20 LH  UPDATES:
% 29 jun 20 LH  changed according to EMGFromLFP

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;

p = inputParser;
addOptional(p, 'basepath', pwd);
addOptional(p, 'fname', '', @ischar);
addOptional(p, 'chunksize', 5e6, @isnumeric);
addOptional(p, 'precision', 'int16', @ischar);
addOptional(p, 'clip', [], @isnumeric);
addOptional(p, 'nchans', 35, @isnumeric);
addOptional(p, 'ch', [33 34 35], @isnumeric);
addOptional(p, 'sensitivity', [0.03573; 0.03551; 0.03455], @isnumeric);
addOptional(p, 'gbias', [0.49445; 0.48464; 0.51740], @isnumeric);
addOptional(p, 'force', false, @islogical);
addOptional(p, 'fsOut', 1250, @isnumeric);
addOptional(p, 'fsIn', 20000, @isnumeric);
addOptional(p, 'graphics', true, @islogical)
addOptional(p, 'saveVar', true, @islogical);

parse(p, varargin{:})
basepath = p.Results.basepath;
fname = p.Results.fname;
chunksize = p.Results.chunksize;
precision = p.Results.precision;
clip = p.Results.clip;
nchans = p.Results.nchans;
ch = p.Results.ch;
sensitivity = p.Results.sensitivity;
gbias = p.Results.gbias;
force = p.Results.force;
fsOut = p.Results.fsOut;
fsIn = p.Results.fsIn;
graphics = p.Results.graphics;
saveVar = p.Results.saveVar;

fprintf('\nextracting accelaration from %s\n', basepath)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% size of one data point in bytes
nbytes = class2bytes(precision);

if isempty(fsOut) % do not resample if new sampling frequency not specified
    fsOut = fsIn;
elseif fsOut ~= fsIn
    [p, q] = rat(fsOut / fsIn);
    n = 5; beta = 20;
    fprintf('accelaration will be resampled from %.1f to %.1f\n', fsIn, fsOut)
end

% handle dat file
cd(basepath)
datfiles = dir([basepath filesep '**' filesep '*dat']);
if isempty(datfiles)
    error('no .dat files found in %s', datpath)
end
if isempty(fname)
    if length(datfiles) == 1
        fname = datfiles.name;
    else
        fname = [bz_BasenameFromBasepath(basepath) '.dat'];
        if ~contains({datfiles.name}, fname)
            error('please specify which dat file to process')
        end
    end
end
[~, basename, ~] = fileparts(fname);

% load acceleration if already exists
accname = [basename '.acceleration.mat'];
if exist(accname, 'file') && ~force
    fprintf('\n loading %s \n', accname)
    load(accname)
    return
end

% partition into chunks
info = dir(fname);
nsamps = info.bytes / nbytes / nchans;
chunks = n2chunks('n', nsamps, 'chunksize', chunksize, 'clip', clip);
nchunks = size(chunks, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load and process acceleration data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% memory map to dat file
m = memmapfile(fname, 'Format', {precision, [nchans, nsamps] 'mapped'});
raw = m.data;

% load timestamps
tname = fullfile(basepath, [basename, '.tstamps.mat']);
if exist(tname, 'file')
    load(tname)
else
    tstamps = 1 : size(m.data.mapped, 2);
end

% go over chunks
a = [];
tstampsACC = [];
for i = 1 : nchunks
    % print progress
    if i ~= 1
        fprintf(repmat('\b', 1, length(txt)))
    end
    txt = sprintf('working on chunk %d / %d', i, nchunks);
    fprintf(txt)
    
    % load data
    d = double(raw.mapped(ch, chunks(i, 1) : chunks(i, 2)));    
    t = double(tstamps(chunks(i, 1) : chunks(i, 2)))';
    
    % convert g to V
    d = (d - gbias) ./ sensitivity;
    
    % remove dc. this is done mainly to diminish edge effects due to the
    % filter embedded within resample (assumes zeros before and after the
    % signal).
    [d] = rmDC(d, 'dim', 2);
    
    % resmaple. tstamps does not require filtering and thus downsample is
    % preferred.
    if fsOut ~= fsIn
        d = [resample(d', p, q, n, beta)]';
        t = downsample(t, q / p);
    end
   
    % calc magnitude (L2 norm)
    mag = vecnorm(d, 2);
    
    % filter
    mag = filterLFP(mag, 'fs', fsOut, 'passband', [0.5 150], 'type', 'butter',...
        'order', 4, 'dataOnly', true, 'graphics', false, 'saveVar', false); 
    
    a = [a; mag];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate power
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% binsize chosen according to FFT in ClusterStates_GetMetrics
binsize = (2 ^ nextpow2(0.5 * fsOut)); 
smf = 7;

% ALT 1: magnitude in 2-5 Hz band (Dhawale et al., eLife, 2017)
% [pband, tband] = specBand('sig', a, 'graphics', false, 'fs', fs,...
%     'band', [2 5], 'binsize', binsize, 'smf', smf, 'normband', false);
% data = bz_NormToRange(pband, [0 1]);

% ALT 2: PC1 of power in broadband
freq = logspace(0, 2, 100);
win = hann(binsize);
[~, ~, tband, pband] = spectrogram(a, win, 0, freq, fsOut, 'yaxis', 'psd');
pband = 10 * log10(abs(pband));
[~, pc1] = pca(pband', 'NumComponents', 1);
data = smooth(pc1, smf);
data = bz_NormToRange(data, [0 1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find episodes of sleep
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fprintf('\nfinding episodes of sleep')
% 
% % find immobility threshold 
% thr = sepBimodel('x', pband, 'lognorm', false, 'nbins', 100,...
%     'smf', 20, 'graphics', false);
% thr = 1400;
% 
% % define episodes of sleep
% vec = [pband < thr];
% sleep = binary2epochs('vec', vec, 'minDur', 10, 'interDur', 1);
% sleep = tband(sleep); % bins to samples
% if sleep == tband(1) % correct edges
%     sleep(1) = 1;
% end
% if ~isempty(sleep)
%     if sleep(end) == tband(end)
%         sleep(end) = tband(length(pband));
%     end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arrange struct and save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
acc.data = data;
acc.tstamps = tband;
acc.fs = fsOut;
acc.fs_orig = fsIn;

if saveVar
    save(accname, 'acc');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if graphics
    figure
    subplot(2, 1, 1)
    plot(acc.tstamps / acc.fs_orig / 60, acc.mag)
    axis tight
    ylabel('L2 Magnitude')
    set(gca, 'TickLength', [0 0], 'XTickLabel', [],...
        'Color', 'none', 'XColor', 'none')
    box off
    subplot(2, 1, 2)
    plot(acc.tband / 60, acc.pband)
    axis tight
    xlabel('Time [m]')
    ylabel('Power Band')
    set(gca, 'TickLength', [0 0], 'Color', 'none')
    box off
    hold on
    Y = ylim;
    fill([acc.sleep fliplr(acc.sleep)]' / 60, [Y(1) Y(1) Y(2) Y(2)],...
        'k', 'FaceAlpha', 0.25,  'EdgeAlpha', 0);
end

fprintf('\nthat took %.2f minutes\n', toc / 60)

end

% EOF