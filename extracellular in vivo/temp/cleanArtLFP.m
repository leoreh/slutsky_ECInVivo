% remove last minutes from recording

basepath = 'E:\Data\Others\DZ\Field\Acute recordings\2h-3h\WT';
cd(basepath)
filename = dir('*.abf');
files = {filename.name};

i = 6;

[~, basename] = fileparts(files{i});

lfp = getLFP('basepath', basepath, 'chans', 1, 'chavg', {},...
    'fs', 1250, 'interval', [0 inf], 'extension', 'abf', 'pli', true,...
    'savevar', true, 'force', false, 'basename', basename);

% remove last x minutes
x = 2;
lfp.data(end : -1 : end - x * 60 * lfp.fs) = [];

% remove first x minutes
x = 2;
lfp.data(1 : fs * 60 * x) = [];

% remove specific artifacts
[~, x] = max(lfp.data);
marg = fs * 0.2;
lfp.data(x - marg : x + marg) = [];

% correct timestamps
lfp.timestamps = [1 : length(lfp.data)] / fs;

% inspect data
figure; plot(lfp.timestamps, lfp.data)


filename = [basename '.lfp.mat'];
save([basepath, filesep, filename], 'lfp')





