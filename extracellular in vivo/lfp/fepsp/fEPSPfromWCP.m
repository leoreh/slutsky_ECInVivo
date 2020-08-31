
%   intens      vec describing stimulus intensity [uA]. must be equal in
%               length to number of recording files in experiment. 
%   dt          numeric. dead time for calculating amplitude. 
%               important for exclusion of stim artifact. {40}[samples]

intens = [100 200 400 800];
stimidx = [40];     % should be in ms according to the protocol, than converted according to fs
fs = 1250;      
inspect = false;
force = true;

%   fs           numeric. requested sampling frequency {1250}



basepath = 'D:\Google Drive\Data\B3\lh61\200826';
cd(basepath)

% load fepsp if already exists
[~, basename, ~] = fileparts(basepath);
fepspname = [basename '.fepsp.mat'];
if exist(fepspname, 'file') && ~force
    fprintf('\n loading %s \n', fepspname)
    load(fepspname)
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% handle files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% find all .wcp files
files = dir('*.wcp');
filenames = natsort({files.name});
% select spicific files
sfiles = [2 : 4];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1 : length(sfiles)
    
    % load data
    filename = files(sfiles(i)).name;
    lfp = getLFP('basepath', basepath, 'basename', filename,...
        'extension', 'wcp', 'forceL', true, 'fs', 1250, 'saveVar', false,...
        'ch', 1, 'cf', 450, 'concat', true, 'dc', true);
    
    % manually inspect and remove unwanted traces
    if inspect
        [lfp.data, rm] = rmTraces(lfp.data, lfp.timestamps);
    end
    
    % arrange data
    wv{i} = lfp.data;
    wvavg(i, :) = mean(lfp.data, 2);
    amp(i) = range(wvavg(i, stimidx : end));
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arrange struct and save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[intens, ia] = sort(intens);

% arrange struct
fepsp.wv = wv(:, ia);
fepsp.wvavg = wvavg(:, ia, :);
fepsp.amp = amp(:, ia);
fepsp.intens = intens;
fepsp.rm = rm;
fepsp.stimidx = stimidx;

if saveVar
    save(fepspname, 'fepsp');
end