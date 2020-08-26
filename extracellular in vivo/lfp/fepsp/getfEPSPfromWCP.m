
%   intens      vec describing stimulus intensity [uA]. must be equal in
%               length to number of recording files in experiment. 
%   dt          numeric. dead time for calculating amplitude. 
%               important for exclusion of stim artifact. {40}[samples]

intens = [100 200 400 800];
dt = 40;

basepath = 'F:\Data\fEPSP\B3\lh59';
cd(basepath)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% handle files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% find all .wcp files
files = dir('*.wcp');
filenames = natsort({files.name});
% select spicific files
sfiles = [1 2 3 5];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1 : length(sfiles)
    
    % load data
    filename = files(sfiles(i)).name;
    lfp = getLFP('basepath', basepath, 'basename', filename,...
        'extension', 'wcp', 'forceL', true, 'fs', 1250, 'saveVar', false,...
        'ch', 1);
    
    % arrange data
    wv{i} = lfp.data;
    wvavg(i, :) = mean(lfp.data, 2);
    amp(i) = range(wvavg(i, dt : end));
    
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
% fepsp.t = win(1) : win(2);
% fepsp.dt = dt;

