function lfp = getLFP(varargin)

% loads lfp data. can specify channels, intervals, average across channels,
% resample, invert. can perform low-pass filter before resampling (see also
% LFPfromDat) and remove DC component. for short repetative signals (e.g.
% stimulation traces) it is best to concatenate signals before filtering
% and removing dc.
%  
% INPUT
%   basename    string. filename of data file. if empty retrieved from
%               basepath. should not include extension.
%   basepath    string. path to load filename and save output {pwd}
%   extension   load from {'lfp'} (neurosuite), 'abf', 'wcp', or 'dat'.
%   forceL      logical. force reload {false}.
%   fs          numeric. requested sampling frequency {1250}
%   interval    numeric mat. list of intervals to read from lfp file [s]
%               can also be an interval of traces from wcp
%   ch          vec. channels to load
%   pli         logical. filter power line interferance {0}
%   dc          logical. remove DC component {0}
%   invertSig   logical. invert signal s.t. max is positive {0}
%   saveVar     logical. save variable {1}.
%   chavg       cell. each row contain the lfp channels you want to average
%   concat      logical. concat traces before removing dc or filtering {1}
%   cf          scalar. cutoff frequency {[450]}.
%   
% DEPENDENCIES
%   import_wcp
%   bz_LoadBinary
%   IOSR.DSP.SINCFILTER
%   LFPfromDat (if extension = 'dat')
%   rmDC
% 
% OUTPUT
%   lfp         structure with the following fields:
%   fs
%   fs_orig
%   extension
%   interval    
%   duration    
%   chans
%   timestamps 
%   data  
% 
% 01 apr 19 LH & RA
% 19 nov 19     load mat if exists  
% 14 jan 19     adapted for wcp and abf 
%               resampling
% 31 aug 20     filter and concat
%
% TO DO LIST
%       # lfp from dat
%       # filter before resampling (done)
%       # load in chunks
%       # downsample by subsampling


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addOptional(p, 'basepath', pwd);
addOptional(p, 'basename', '');
addOptional(p, 'extension', 'lfp');
addOptional(p, 'forceL', false, @islogical);
addOptional(p, 'fs', 1250, @isnumeric);
addOptional(p, 'interval', [0 inf], @isnumeric);
addOptional(p, 'ch', [1 : 16], @isnumeric);
addOptional(p, 'pli', false, @islogical);
addOptional(p, 'dc', false, @islogical);
addOptional(p, 'invertSig', false, @islogical);
addOptional(p, 'saveVar', true, @islogical);
addOptional(p, 'chavg', {}, @iscell);
addOptional(p, 'concat', false, @islogical);
addOptional(p, 'cf', 450, @isnumeric);

parse(p,varargin{:})
basepath = p.Results.basepath;
basename = p.Results.basename;
extension = p.Results.extension;
forceL = p.Results.forceL;
fs = p.Results.fs;
interval = p.Results.interval;
ch = p.Results.ch;
pli = p.Results.pli;
dc = p.Results.dc;
invertSig = p.Results.invertSig;
saveVar = p.Results.saveVar;
chavg = p.Results.chavg;
concat = p.Results.concat;
cf = p.Results.cf;

nchans = length(ch);

if isempty(basename)
    [~, basename] = fileparts(basepath);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check if file exists
cd(basepath)
filename = [basename '.lfp.mat'];
if exist(filename) && ~forceL
    fprintf('\n loading %s \n', filename)
    load(filename)
    return
end

loadname = [basename '.' extension];
fprintf('\nworking on %s\n', loadname)
switch extension
    case 'lfp'
        fs_orig = 1250;
        sig = bz_LoadBinary(loadname, 'duration', diff(interval),...
            'frequency', fs_orig, 'nchannels', nchans, 'start', interval(1),...
            'channels', ch, 'downsample', 1);
    case 'abf'
        % note abf2load cannot handles spaces in loadname
        % note abf2load requires Abf2Tmp.exe and ABFFIO.dll in basepath         
        [sig, info] = abf2load(loadname);
        fs_orig = 1 / (info.fADCSequenceInterval / 1000000); 
    case 'wcp'
        data = import_wcp(loadname);
        data.S = [data.S{ch}];
        if interval(2) ~= Inf
            data.S = data.S(:, interval(1) : interval(2));
        end
        sig = data.S;
        fs_orig = data.fs;
    case 'dat'
        error('not ready yet')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% messaround
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% low-pass filter
if cf
    fprintf('low-pass filtering, cutoff = %d Hz\n', cf)
    import iosr.dsp.*
    filtRatio = cf / (fs_orig / 2);
    if concat
        sz = size(sig);
        sig = sig(:);
    end
    sig = iosr.dsp.sincFilter(sig, filtRatio);
    if concat
        lfp.data = reshape(sig, sz);
    end
end

% remove DC component
if dc
    fprintf('removing dc component\n')
    lfp.data = rmDC(lfp.data, 'dim', 1);
end

% resmaple
if isempty(fs) % do not resample if new sampling frequency not specified 
    fs = fs_orig;
elseif round(fs) ~= round(fs_orig)
    fprintf('resampling from %.1f to %.1f Hz\n', fs_orig, fs)
    [p, q] = rat(fs / fs_orig);
    n = 5; beta = 20;
    for i = 1 : size(lfp.data, 2)        % only way to handle large arrays
        draw(:, i) = [resample(double(lfp.data(:, i))', p, q, n, beta)]';
    end
    lfp.data = draw;    
end
lfp.timestamps = (interval(1) : 1 / fs : interval(1) + (length(lfp.data) - 1) / fs)';

% invert
if invertSig
    if ~strcmp(extension, 'lfp') && abs(min(lfp.data(:))) > max(lfp.data(:))
        fprintf('\n inverting data \n\n')
        lfp.data = -lfp.data;
    end
end

% convert to double
lfp.data = double(lfp.data);

% filter power line interference
if pli
    linet = lineDetect('x', lfp.data, 'fs', fs, 'graphics', false);
    lfp.data = lineRemove(lfp.data, linet, [], [], 0, 1);
end

% flip such that samples x channels
if size(lfp.data, 1) < size(lfp.data, 2)
    lfp.data = lfp.data';
end

% signal average
if ~isempty(chavg)
    mlfp = zeros(size(chavg, 1), length(lfp.data));
    for i = 1 : size(chavg, 1)
        mlfp(i, :) = mean(lfp.data(:, chavg{i}), 2);
    end    
    lfp.data = mlfp;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arrange and save struct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if interval(2) == inf
    inteval(2) = lfp.timestamps(end);
end
lfp.interval = interval;
lfp.duration = length(lfp.data) / fs;
lfp.chans = ch;
lfp.fs = fs;
lfp.fs_orig = fs_orig;
lfp.origFile = extension;

% save variable
if saveVar   
    save([basepath, filesep, filename], 'lfp')
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sanity checks and testing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for trace = 1 : size(lfp.data, 2)
    figure
    plot(lfp.timestamps, lfp.data(:, trace));
    hold on
    plot(data.T, data.S(:, trace))  
end

% conclusion from below: no difference if downsample average or average the
% downsampled signal
figure
plot(lfp.timestamps, mean(lfp.data, 2))
hold on
plot(data.T, mean(data.S, 2))
dsamp = resample(mean(data.S, 2), p, q, n, beta);
plot(lfp.timestamps, dsamp)
legend

% filter
% cf = 450;
% filtRatio = cf / (fs_orig / 2);
% filtered = [iosr.dsp.sincFilter(data.S, filtRatio)];
% 
% for trace = 2 : size(lfp.data, 2)
%     figure
%     plot(data.T, filtered(:, trace))
%     hold on
%     plot(data.T, data.S(:, trace))
% end

% resample. both alternatives below produce similar results, alt 2 is
% faster
% ALT 1: matlab resmpale with antialiasing filter
[p, q] = rat(fs / lfp.fs_orig);
n = 5; beta = 20;
sig = resample(sig, p, q, n, beta);

% ALT 2: downsample by subsampling
fsRatio = lfp.fs_orig / fs;
sig = sig(1 : fsRatio : length(sig), :);
t = 0 : 1 / fs : ((length(sig2) - 1) / fs);