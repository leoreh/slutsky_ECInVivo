
% traces params
stimN = 5;
ch = 3;
fs = 20000;
intens = [20 : 20 : 100];

% file params
protocol_id = 'io';
basepath = 'I:\Data\Processed\lh50\fepsp\lh50_200426_120451';
[~, basename] = fileparts(basepath);
cd(basepath)

% load data
data_in = double(bz_LoadBinary([basename, '.dat'], 'frequency', fs, 'nChannels', 31,...
    'channels', [11, 13, 21]));

% convert tstamps to idx of samples
load([basename, '.din.mat'])
load([basename, '.tstamps.mat'])
stimVec = zeros(1, length(din.data));
for istim = 1 : length(din.data)
    stimVec(istim) = find(tstamps == din.data(istim));
end
        
% check
fh = figure;
plot(data_in(:, ch))
hold on
plot([stimVec,; stimVec], ylim, '--k')


% org stimIdx to cell of intens. use file duration
load([basename, '.datInfo.mat'])
nfiles = length(datInfo.origFile);
nsamps = datInfo.nsamps;
csamps = [0 cumsum(nsamps)];
maxstim = 1;
stim_locs = cell(nfiles, 1);
for ifile = 1 : nfiles
    stim_locs{ifile} = stimVec((stimVec > csamps(ifile) &...
        stimVec <  csamps(ifile + 1)));
end

% save to new file
x(1, :) = -data_in(:, 1);
x(2, :) = -data_in(:, 2);
x(3, :) = data_in(:, 3);
x = int16(x);
basename = 'fepsp_testData';
dataname = [basename, '.dat'];
fid = fopen(dataname, 'w');
fwrite(fid, x(:), 'int16');
fclose('all')

save([basename, '_stimLocs.mat'], 'stim_locs')