function data = binary_load(filename, varargin)
% Reads data from a binary file, allowing for selection of specific channels,
% time ranges or sample ranges, and downsampling.
%
% INPUT
%   filename    - char array, path to the binary file.
%   fs          - numeric, sampling rate in Hz (default: 20000).
%   start       - numeric, position to start reading in seconds (default: 0).
%                 Cannot be used with offset or samples.
%   duration    - numeric, duration to read in seconds (default: Inf, reads to end).
%                 Cannot be used with offset or samples.
%   offset      - numeric, position to start reading in samples per channel
%                 (default: 0). Cannot be used with start or duration.
%   samples     - numeric, number of samples per channel to read 
%                 (default: Inf, reads to end). Cannot be used with start or duration.
%   nCh         - numeric, total number of channels in the file (default: 1).
%   ch          - numeric vector, specific channels to read (1-based, 
%                 default: all channels from 1 to nCh).
%   precision   - char array, sample precision (e.g., 'int16', 'double', 
%                 default: 'int16'). See `fread` documentation for valid types.
%   skip        - numeric, number of bytes to skip after each frame of 
%                 nCh samples is read (default: 0).
%                 Note: This skip is applied *after* a full block of 
%                 `nCh` samples conforming to `precision` is read.
%                 If `downsample > 1`, this parameter is ignored, and skipping
%                 is determined by the downsample factor.
%   downsample  - numeric, factor by which to downsample (default: 1, no downsampling).
%                 Must be an integer >= 1.
%   bit2uv      - numeric scalar or empty array. Conversion factor to microvolts.
%                 Behavior is as follows:
%                   - If 'bit2uv' is explicitly passed as empty (`[]`), `bit2uv` is set to 1.
%                   - Else if `round(fs) == 24414` (e.g., common TDT sampling rate), `bit2uv` is set to 1.
%                   - Else if 'bit2uv' is explicitly passed as a non-empty scalar, that value is used.
%                   - Else (if 'bit2uv' is not provided by the user), it defaults to 0.195.
%
% OUTPUT
%   data        - matrix (nSamples x numel(ch)) of loaded data.
%                 The data type is `double`.
%                 The data is scaled by the `bit2uv` factor *twice*.
%
% Based on bz_LoadBinary by Michaël Zugaro (2004-2011) and DLevenstein (2016).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parse Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;

% Required
addRequired(p, 'filename', @(x) ischar(x) || isstring(x));

% Optional Parameters
addParameter(p, 'fs', 20000, @(x) isnumeric(x) && isscalar(x) && isreal(x) && x > 0);
addParameter(p, 'start', 0, @(x) isnumeric(x) && isscalar(x) && isreal(x) && x >= 0);
addParameter(p, 'duration', Inf, @(x) isnumeric(x) && isscalar(x) && isreal(x) && x >= 0);
addParameter(p, 'offset', 0, @(x) isnumeric(x) && isscalar(x) && isreal(x) && mod(x,1)==0 && x >= 0);
addParameter(p, 'samples', Inf, @(x) isnumeric(x) && isscalar(x) && isreal(x) && x >= 0); 
addParameter(p, 'nCh', 1, @(x) isnumeric(x) && isscalar(x) && isreal(x) && mod(x,1)==0 && x > 0);
addParameter(p, 'ch', [], @(x) (isnumeric(x) && isvector(x) && isreal(x) && all(mod(x,1)==0) && all(x > 0)) || isempty(x));
addParameter(p, 'precision', 'int16', @(x) ischar(x) || isstring(x));
addParameter(p, 'skip', 0, @(x) isnumeric(x) && isscalar(x) && isreal(x) && mod(x,1)==0 && x >= 0);
addParameter(p, 'downsample', 1, @(x) isnumeric(x) && isscalar(x) && isreal(x) && mod(x,1)==0 && x >= 1);
addParameter(p, 'bit2uv', 0.195, @(x) isnumeric(x) && (isscalar(x) || isempty(x)));

parse(p, filename, varargin{:});

% Assign parsed inputs to local variables
filename = char(p.Results.filename); 
fs = p.Results.fs;
start = p.Results.start;
duration = p.Results.duration;
offset = p.Results.offset;
samples = p.Results.samples;
nCh = p.Results.nCh;
ch = p.Results.ch;
precision = char(p.Results.precision);
skip = p.Results.skip;
downsample = p.Results.downsample;
bit2uv = p.Results.bit2uv;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Handle bit2uv default logic (for tdt)
if isempty(bit2uv) || round(fs) == 24414
    bit2uv = 1;
end

% Handle ch default
if isempty(ch)
    ch = 1:nCh;
end

% Check for parameter conflicts
params_passed = setdiff(p.Parameters, p.UsingDefaults);
timeParams = any(ismember({'start', 'duration'}, params_passed));
sampleParams = any(ismember({'offset', 'samples'}, params_passed));

if timeParams && sampleParams
    error('Data params specified both in time and in samples.');
end

% Check consistency between channel IDs and number of channels
if any(ch > nCh) || any(ch < 1)
    error('Cannot load specified channels.');
end

% Determine sampleSize in bytes
sampleSize = class2bytes(precision);
if isempty(sampleSize)
    error('Unsupported precision string: %s.', precision);
end

% Calculate dataOffset_bytes and nSamples based on specification mode
if timeParams
    offsetBytes = floor(start * fs) * nCh * sampleSize;
    nSamples = round(duration * fs);
else
    offsetBytes = offset * nCh * sampleSize;
    nSamples = samples;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File Operations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Open file
if ~exist(filename, 'file')
    error('File ''%s'' not found.', filename);
end
f = fopen(filename,'r');
if f == -1
    error('Cannot read ''%s''.', filename);
end

% Position file index for reading
status = fseek(f,offsetBytes,'bof');
if status ~= 0
    fclose(f);
    error('Could not start reading at specified offset (%.0f bytes).', offsetBytes);
end

% Determine total number of samples available in file from current position
fileCurrentPos = ftell(f);
fseek(f,0,'eof'); 
fileStopPos = ftell(f);
fseek(f,fileCurrentPos,'bof'); 

maxBytes = fileStopPos - fileCurrentPos;
if maxBytes < 0
    maxBytes = 0;
end 
maxSamples = floor(maxBytes / (nCh * sampleSize));

if isinf(nSamples) || nSamples > maxSamples
    nSamples = maxSamples;
end

% Handle downsampling and skip parameter for LoadChunk
if downsample > 1
    skipSamples = nCh * (downsample - 1) * sampleSize;
    nSamples = floor(nSamples / downsample);
else
    skipSamples = skip; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Determine chunks for reading using nSamples (total samples after downsampling/limits)
chunks = n2chunks('n', nSamples, 'chunksize', 100000); 
nChunks = size(chunks, 1);

if nChunks == 1
    chunkSamples = chunks(1, 2) - chunks(1, 1) + 1;
    data = LoadChunk(f, nCh, ch, chunkSamples, precision, skipSamples);
    
    % Apply bit2uv conversion
    data = double(data) * bit2uv; 
else 
    
    % Preallocate data matrix with the file's precision.
    data = zeros(nSamples, numel(ch), precision);
    samplesCnt = 0;

    for iChunk = 1:nChunks
        chunkSamples = chunks(iChunk, 2) - chunks(iChunk, 1) + 1;   

        chunkData = LoadChunk(f, nCh, ch, chunkSamples, precision, skipSamples);
   
        idxStart = samplesCnt + 1;
        idxEnd = samplesCnt + chunkSamples;
        
        data(idxStart : idxEnd, :) = double(chunkData) * bit2uv;
        samplesCnt = idxEnd;
        
        % The condition 'chunkSamples < chunkSamples' is always false.
        % Thus, this break statement, intended for EOF detection within a chunk, is currently not triggered.
        if chunkSamples < chunkSamples
            break; 
        end
    end
end

% Assure all data read
if samplesCnt < nSamples
   error('Less data than expected read')
end

% close file
fclose(f);

% Apply bit2uv conversion after all data is loaded and matrix is finalized
    data = double(data) * bit2uv; 

end % EOF

% ---------------------------------------------------------------------------------------------------------
% Local Helper Function
% ---------------------------------------------------------------------------------------------------------
function data = LoadChunk(fid, nCh, ch, chunkSamples, precision, skipSamples)
% LoadChunk - Reads a chunk of data from the binary file.
%
% INPUTS:
%   fid          - File identifier for the open binary file.
%   nCh          - Total number of channels multiplexed in the file.
%   ch           - Vector of 1-based channel indices to extract.
%   chunkSamples - Number of samples per channel to read in this chunk.
%   precision    - String defining the data precision (e.g., 'int16').
%   skipSamples  - Number of bytes to skip in the file after reading each
%                  frame of `nCh` samples. This is used either for user-defined
%                  skipping or for implementing downsampling.
%
% OUTPUT:
%   data         - Matrix (actual_samples_read_in_chunk x numel(ch))
%                  containing the loaded data for the specified channels for this chunk.
%                  `actual_samples_read_in_chunk` might be less than `chunkSamples` if EOF is reached.

% Format string for fread: read a full frame of all channels
readFormat = sprintf('%d*%s=>%s', nCh, precision, precision);

if skipSamples ~= 0
    data = fread(fid, [nCh, chunkSamples], readFormat, skipSamples);
else
    data = fread(fid, [nCh, chunkSamples], readFormat);
end

% Transpose to have samples as rows, channels as columns
data = data';

% Extract specified channels
data = data(:, ch);

end 