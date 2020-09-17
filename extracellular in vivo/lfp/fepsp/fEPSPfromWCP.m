function fepsp = fEPSPfromWCP(varargin)

% wrapper to get fEPSP signals recorded via winWCP. assumes folder
% (basepath) contains several homogenous files (i.e. each file consists of
% only one type of stimulus). files are analyzed according to stimulus
% protocol (e.g. io or stp). user is asked to input / select specific
% files. data is loaded via getLFP such that traces are concatenated,
% low-pass filtered and downsampled. 
%  
% INPUT
%   basepath    string. path to data files. if sfiles is specified than
%               basepath needs be to to the folder with wcp files
%   sufx        string to add to filename (e.g. 'io1') {''}
%   intens      vec describing stimulus intensity [uA]. must be in the same
%               order as the recording files
%   protocol    stimulus protocol. can be a string ('io' or 'stp') or a
%               numeric vector describing the times [ms] of stimulations
%               for each trace (currently not implemented)
%   dt          numeric. dead time after stimulus (for exluding stim
%               artifact) {2}[ms]
%   sfiles      vec of file indices to analyse. if empty a prompt will ask
%               user to select files
%   fs          numeric. requested sampling frequency {1250}
%   inspect     logical. inspect traces {0}.
%   saveVar     logical. save variable {1}.
%   saveFig     logical. save graphics {1}. 
%   graphics    logical. plot graphics {1}. 
%   force       logical. force reload {0}.
% 
% DEPENDENCIES
%   getLFP
%   cell2nanmat
%   rmTraces
%   
% OUTPUT
%   fepsp       struct
% 
% TO DO LIST
%   # implement arbitrary stim protocols
%   # analyse more parameters than amp
% 
% 02 sep 20 LH   UPDATES:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addOptional(p, 'basepath', pwd);
addOptional(p, 'sufx', '', @ischar);
addOptional(p, 'sfiles', [], @isnumeric);
addOptional(p, 'intens', [], @isnumeric);
addOptional(p, 'protocol', 'io', @ischar);
addOptional(p, 'dt', 2, @isnumeric);
addOptional(p, 'fs', 1250, @isnumeric);
addOptional(p, 'force', false, @islogical);
addOptional(p, 'saveVar', true, @islogical);
addOptional(p, 'saveFig', true, @islogical);
addOptional(p, 'graphics', true, @islogical);
addOptional(p, 'inspect', false, @islogical);

parse(p, varargin{:})
basepath = p.Results.basepath;
sufx = p.Results.sufx;
sfiles = p.Results.sfiles;
intens = p.Results.intens;
protocol = p.Results.protocol;
dt = p.Results.dt;
fs = p.Results.fs;
force = p.Results.force;
saveVar = p.Results.saveVar;
saveFig = p.Results.saveFig;
graphics = p.Results.graphics;
inspect = p.Results.inspect;

ylimit = [-1.5 1.5];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt = ceil(dt / 1000 * fs);      % ms to samples

cd(basepath)
% find all .wcp files
if ~isempty(sfiles)
    files = dir('*.wcp');
    filenames = natsort({files.name});
    filenames = filenames(sfiles);
else
    [filenames, basepath] = uigetfile('*.wcp', 'MultiSelect', 'on');
end

% session metadata; assumes file structure: animal/date/
pathPieces = regexp(basepath, filesep, 'split');
if isempty(pathPieces{end})
    pathPieces(end) = [];
end
mouse = pathPieces{end - 1};
expDate = pathPieces{end};
expName = [mouse '_' expDate];
if ~isempty(sufx)
    expName = [expName '_' sufx];
end

% load fepsp if already exists
fepspname = [expName '.fepsp.mat'];
if exist(fepspname, 'file') && ~force
    fprintf('\n loading %s \n', fepspname)
    load(fepspname)
    return
end

% get timing of stimulus in samples
switch protocol
    case 'io'
        % single pulse of 500 us after 30 ms. recording length 150 ms.
        % repeated once every 15 s
        nstim = 1;
        stimidx = ceil(30 / (1 / fs * 1000));
    case 'stp'
        % 5 pulses of 500 us at 50 Hz. starts after 10 ms. recording length
        % 200 ms. repeated once every 30 s
        tstim = 1 / 50 * 1000;      % ms
        nstim = 5;
        stimidx = [10 : tstim : (nstim - 1) * tstim + 10]; 
        stimidx  = ceil(stimidx / (1 / fs * 1000));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1 : length(filenames)
    
    % load data
    [~, basename] = fileparts(filenames{i});
    lfp = getLFP('basepath', basepath, 'basename', basename,...
        'extension', 'wcp', 'forceL', true, 'fs', fs, 'saveVar', false,...
        'ch', 1, 'cf', 450, 'concat', true, 'dc', true);
    
    lfp.data = -lfp.data;
    
    % manually inspect and remove unwanted traces
    if inspect
        [lfp.data, rm] = rmTraces(lfp.data, 'x', lfp.timestamps,...
            'ylim', ylimit);
        lfp.data(:, rm) = [];
    else
        rm = [];
    end
    
    % analyze data; this is your part lior :)
    wv{i} = lfp.data;
    wvavg(i, :) = mean(lfp.data, 2);
    
    switch protocol
        case 'io'
            amp(i) = range(wvavg(i, stimidx + dt : end));
            ampcell{i} = range(wv{i}(stimidx + dt : end, :));
        case 'stp'
            for ii = 1 : nstim
                s1 = stimidx(ii) + dt;
                s2 = floor(min([stimidx(ii) + tstim / 1000 * fs - dt, length(wv{i})]));
                amp(i, ii) = range(wvavg(i, s1 : s2));
                ampcell{i}(ii, :) = range(wv{i}(s1 : s2, :));
            end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arrange struct and save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[intens, ia] = sort(intens);

% arrange struct
fepsp.origFiles = filenames;
fepsp.wv = wv(:, ia);
fepsp.timestamps = lfp.timestamps;
fepsp.wvavg = wvavg(ia, :);
fepsp.amp = amp(ia);
fepsp.ampmat = cell2nanmat(ampcell(ia));
fepsp.intens = intens;
fepsp.rm = rm;
fepsp.stimidx = stimidx;
fepsp.fs = lfp.fs;

if saveVar
    save(fepspname, 'fepsp');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if graphics
    figure
    suptitle(expName)
    subplot(1, 2, 1)    % waveform
    hold on
    for i = 1 : length(intens)
        plot(fepsp.timestamps, fepsp.wvavg(i, :));
    end
    legend(split(num2str(intens)))
    xlabel('Time [ms]')
    ylabel('Voltage [mV]')
    box off
    ylim(ylimit)

    subplot(1, 2, 2)    % amplitude
    switch protocol
        case 'io'
            boxplot(fepsp.ampmat, 'PlotStyle', 'traditional')
            xticklabels(split(num2str(intens)))
            xlabel('Intensity [uA]')
        case 'stp'
            for i = 1 : length(filenames)
                plot(amp(i, :) / amp(i, 1))
                hold on
                ylim([0 1])
                xlabel('Stim. No.')
            end
    end
    ylabel('Amplidute [mV]')
    box off
    
    if saveFig
        export_fig(expName, '-tif', '-r300', '-transparent')
    end
end