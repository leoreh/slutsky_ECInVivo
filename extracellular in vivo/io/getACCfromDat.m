function din = getACCfromDat(varargin)

% get 3-axis accelerometer signal from dat file and return the acceleration
% amplitude. based on
% http://intantech.com/files/Intan_RHD2000_accelerometer_calibration.pdf
% and Dhawale et al., eLife, 2017. Conversion from g to volts is done by
% substracting the zero-g bias [V] and dividing with the accelerometer
% sensitivity [V / g]. these parameters are based on measurments done
% XX.XX.XX. After conversion to Volts, the 3-axis acceleration data is L2
% normalized to yeild a magnitude heuristic. 
%
% INPUT:
%   basepath    string. path to .dat and .npy files {pwd}. 
%               if multiple dat files exist than concat must be true
%   fname       string. name of dat file. if empty and more than one dat in
%               path, will be extracted from basepath
%   precision   char. sample precision {'int16'} 
%   nchans      numeric. number of channels in dat file {35}.
%   ch          3x1 vec. channel index to acceleration in x,
%               y, z (order does not matter so long corresponds to
%               sensitivity and gbias).
%   sensitivity 3x1 vec. sensitivity data to acceleration in x,
%               y. order must correspond to ch. 
%   gbias       3x1 vec. zero-g bias data to acceleration in x,
%               y. order must correspond to ch.  
%   force       logical. reload data even if .mat file exists {false}
%
% OUTPUT
%   din         
%
% CALLS:
%   class2bytes
%
% TO DO LIST:
%   # downsampling
%   # linear envelop (smoothing and filtering)
%
% 09 apr 20 LH      

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addOptional(p, 'basepath', pwd);
addOptional(p, 'fname', '', @ischar);
addOptional(p, 'precision', 'int16', @ischar);
addOptional(p, 'nchans', 35, @isnumeric);
addOptional(p, 'ch', [33 34 35], @isnumeric);
addOptional(p, 'sensitivity', [0.3468; 0.3468; 0.3468], @isnumeric);
addOptional(p, 'gbias', [1.775; 1.775; 1.775], @isnumeric);
addOptional(p, 'force', false, @islogical);

parse(p, varargin{:})
basepath = p.Results.basepath;
fname = p.Results.fname;
precision = p.Results.precision;
nchans = p.Results.nchans;
ch = p.Results.ch;
sensitivity = p.Results.sensitivity;
gbias = p.Results.gbias;
force = p.Results.force;

basepath = 'E:\Data\Dat\lh50\lh50_200402\190448_e2r1-7';
precision = 'int16';
nchans = 31;
nbytes = class2bytes(precision);
ch = [29 30 31];
sensitivity = [0.3468; 0.3468; 0.3468]; % [V/g]
gbias = [1.775; 1.775; 1.775]; % V

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
accname = [basename '.acc.mat'];
if exist(accname, 'file') && ~force
    fprintf('\n loading %s \n', accname)
    load(accname)
    return
end

% memory map to dat file
info = dir(fname);
nsamps = info.bytes / nbytes / nchans;
m = memmapfile(fname, 'Format', {precision, [nchans, nsamps] 'mapped'});

% load data
acc = double(m.Data.mapped(ch, :));

% resmaple
if isempty(fs) % do not resample if new sampling frequency not specified 
    fs = fs_orig;
elseif fs ~= fs_orig
    [p, q] = rat(fs / fs_orig);
    n = 5; beta = 20;
    for i = 1 : size(lfp.data, 2)        % only way to handle large arrays
        draw(:, i) = [resample(double(lfp.data(:, i))', p, q, n, beta)]';
    end
    lfp.data = draw;    
    fprintf('\n resampling from %.1f to %.1f\n\n', fs_orig, fs)
end

% convert g to V
acc = (acc - gbias) ./ sensitivity;

% downsample

% calc magnitude (L2 norm)
mag = vecnorm(acc, 2);

% filter
mag = filterLFP(mag, 20000, 'passband', [0.5 150], 'type', 'butter',...
    'order', 4, 'dataOnly', true, 'graphics', false, 'saveVar', false);

% arrange struct output and save
din.data = data;
din.origDir = {jsonFiles.folder};
if saveVar
    save(destination, 'din');
end

end

% EOF