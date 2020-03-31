
% this is a wrapper to inspect and clean lfp signals

% user argument
basepath = 'E:\Data\Field\IIS\APPPS1';
i = 5;

% file
[~, grpname] = fileparts(basepath);
cd(basepath)
filename = dir('*.lfp.*');
files = natsort({filename.name});
[~, basename] = fileparts(files{nfiles(i)});
[~, basename] = fileparts(basename);
    
% parameters    
binsize = (2 ^ nextpow2(30 * 1250));
marg = round(fs * 0.05);        % rm 50 ms arround an artifact

% load lfp
lfp = getLFP('basepath', basepath, 'ch', 1, 'chavg', {},...
    'fs', 1250, 'interval', [0 inf], 'extension', 'abf', 'pli', true,...
    'dc', true, 'savevar', true, 'force', false, 'basename', basename);
                   
% find and remove bad episodes
ep = markEp(lfp.timestamps / 60, lfp.data);
for i = 1 : size(ep, 1)
    lfp.data(ep(i, 1) : ep(i, 2) = [];
end

% inspect data
figure; plot(lfp.timestamps / 60, lfp.data)

% invert 
lfp.data = -lfp.data;

% remove max artifacts
[~, x] = max(lfp.data);
lfp.data(x - marg : x + marg) = [];

% correct timestamps
lfp.timestamps = [1 : length(lfp.data)] / lfp.fs;

% save
filename = [basename '.lfp.mat'];
save([basepath{i}, filesep, filename], 'lfp')
    


