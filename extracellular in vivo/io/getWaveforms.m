function wav = getWaveforms(filename, varargin)

% get raw waveform during specfic times
% 
% INPUT:
%   filename    full name of spk file including path
%   clu         clusters of interest. given as a cluID vector. if empty
%               all clusters will be loaded
%   spkidx      spikes of interest. given as a cell of idx vectors. the
%               number of cells must be equal to length(clu). for
%               example, to get the waveform of spikes fired after 5
%               min; spkidx = find(spikes.times > 300). spkidx must be
%               a cell even if only one cluster if of interest
%   nsamps      number of samples per waveform {32}.
%   nchans      number of channels in spike group {4}.
%   saveVar     save output in basepath
%
% OUTPUT:
%   wav         array of k cells where k is the number of clusters. each
%               cell is a matrix of n x m x l, where n is the number of
%               spikes, m is the number of electrodes, and l is the number
%               of samples.
%
% 16 jan 19 LH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addOptional(p, 'clu', []);
addOptional(p, 'spkidx', [], @iscell);
addOptional(p, 'nsamps', 32, @isnumeric);
addOptional(p, 'nchans', 4, @isnumeric);
addOptional(p, 'saveVar', false, @islogical);

parse(p,varargin{:})
clu = p.Results.clu;
spkidx = p.Results.spkidx;
nsamps = p.Results.nsamps;
nchans = p.Results.nchans;
saveVar = p.Results.saveVar;

if length(clu) ~= length(spkidx)
    error('number of cells in spkidx must be equal to number of clusters')
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[basepath, ~] = fileparts(filename);
cd(basepath)

% open spk file and load waveforms
fid = fopen(filename, 'r');
if(fid == -1)
    error('cannot open file');
end
rawwav = fread(fid, [1 inf], 'int16=>int16');
rawwav = reshape(rawwav, nchans, nsamps, []);
rawwav = permute(rawwav,[3 1 2]);
fclose(fid);

% open corresponding clu file
[basepath, ~, grp] = fileparts(filename);
[~, basename] = fileparts(basepath);
cluname = [basename, '.clu', grp];
fid = fopen(cluname, 'r');
if(fid == -1)
    error('cannot open file');
end
nclu = fscanf(fid, '%d', 1);
cluvec = fscanf(fid, '%f')';
fclose(fid);
if isempty(clu)
    clu = unique(cluvec);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% select data of interest
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1 : length(clu)
    % select clusters
    wav{i} = rawwav(cluvec == clu(i), :, :);
    % select spikes
    wav{i} = wav{i}(spkidx{i}, :, :);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if saveVar
    save(fullfile(basepath, basename, 'wav.mat'), 'wav')
end

end

% EOF
