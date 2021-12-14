function fepsp = fEPSPfromWCP(varargin)

% wrapper to get fEPSP signals recorded via winWCP. assumes folder
% (basepath) contains several homogenous files (i.e. each file consists of
% only one type of stimulus). files are analyzed according to stimulus
% protocol (e.g. io or stp). user is asked to input / select specific
% files. data is loaded via getLFP such that traces may be concatenated,
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
%   dt          numeric. dead time before / after stimulus (for exluding stim
%               artifact) {2}[ms]
%   sfiles      vec of file indices to analyse. if empty a prompt will ask
%               user to select files
%   fs          numeric. requested sampling frequency {1250}
%   inspect     logical. inspect traces {0}. obsolete after Lilu's gui
%   saveVar     logical. save variable {1}.
%   saveFig     logical. save graphics {1}. Only relevant If anaflag is
%               true
%   graphics    REMOVED! - logical. plot graphics {1}.
%   force       logical. force reload {0}.
%   anaflag     logical. send to analysis {1}
%   AddAnalyseParm
%               Additional Parameters to pass to fEPSP_analysis, as cell
%               vector (see fEPSP_analysis for more info). { {} }
% 
% DEPENDENCIES
%   getLFP
%   cell2nanmat
%   rmTraces
%   fEPSP_analysis
%   
% OUTPUT
%   fepsp           struct with fields:
%       Info        struct, info about the recourding, with fields:
%           basename    cell, the name of the file/ files that the traces were
%                       exstracted from.
%           recSystem   char, type of recourding system (wcp/oe/tdt), in
%                       this function always 'wcp'.
%           protocol    stimulus protocol. can be a string ('io' or 'stp') or a
%                       numeric vector describing the times [ms] of stimulations
%                       for each trace (numeric currently not implemented)
%           fsOrig      double, original sampling frequency of file.
%           fs          double, sampling frequency after downsampling
%           spkgrp      cell, for which channles are realated to which
%                       electrode. For this function always 1 electrode
%                       countaining 1 channle, so {1}.
%           stimIdx     double, the indexes of the stimuli in the long dat
%                       file (in this function empty and for compatability
%                       with fEPSP from dat only) (need to double check meaning).
%           stimTs      double, the index of the stimuli in the TS of the
%                       long dat file (in this function empty and for
%                       compatability with fEPSP from dat only).
%           stamps      double, the index of the start and end of each
%                       trace in the long dat file (in this function empty
%                       and for compatability with fEPSP from dat only).
%           intensOrig  double, the intensities as user inputed them, without 
%                       any sorting and merging.
%           rm          cell, for each intensity the number of the traces
%                       user removed. Column order matching intensOrig.
%           lowPass     double, the cut off frequency for the low pass filter used.
%           inspect     logical, true if user inspected the traces in each
%                       intensity in order to remove unwanted ones.
%       intens      double the intensities in this file after sorting only
%                   in this function, and also merging of duplicated in fEPSPfromDat.
%       tstamps     the time stamps for all traces. 0 is the first stimulus.
%       traces      cell, all the traces extracted, electrode X intensities
%                   (electrode only 1 in this function).
% TO DO LIST
%   # implement arbitrary stim protocols
%   # analyse additional parameters (slope, pSpike) # Moved to fEPSP_analysis
% 
% 02 sep 20 LH   UPDATES:
% 26 sep 20 LH      adaptations according to clampfit
% 31 oct 20 LD      match output epsp format to the output from
%                       fEPSPfromDat, and transfer analyse to fEPSP_analysis
% 05 dec 20 LD      When anaflag is true, fepsp output + saved will be the
%                   output of fEPSP_analysis

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
%addOptional(p, 'graphics', true, @islogical);
addOptional(p, 'inspect', false, @islogical);
addOptional(p, 'anaflag', true, @islogical);
addOptional(p, 'AddAnalyseParm',{},@(x) validateattributes(x,"cell",{'vector'}))

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
%graphics = p.Results.graphics;
inspect = p.Results.inspect;
anaflag = p.Results.anaflag;
AddAnalyseParm = p.Results.AddAnalyseParm;  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt = ceil(dt / 1000 * fs);      % ms to samples
ylimit = [-1.5 1.5];

cd(basepath)
% find all .wcp files
if ~isempty(sfiles)
    files = dir('*.wcp');
    filenames = natsort({files.name});
    filenames = filenames(sfiles);
else
    [filenames, basepath] = uigetfile('*.wcp', 'MultiSelect', 'on');
    if ~iscell(filenames)
        filenames = {filenames};
    elseif isempty(filenames)
        error('No files selected. Aborting')
    end
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
    if anaflag
        fepsp = fEPSP_analysis('fepsp', fepsp,'saveFig',saveFig,'saveVar',saveVar,'savename',fepspname);
    end
    return
end

% get timing of stimulus in samples
switch protocol
    case 'io'
        % single pulse of 500 us after 30 ms. recording length 150 ms.
        % repeated once every 15 s. negative peak of response typically
        % 10 ms after stim.
%         nstim = 1;
%         stimidx = ceil(30 / (1 / fs * 1000)) + dt;        % sample of stim exlcuding deadtime
        baseline = [1 floor(30 * fs / 1000) - dt];        % samples of baseline period
%         ampidx = [stimidx : stimidx + 15 * fs / 1000];    % samples for calculting amplitude
%         wvidx = [stimidx : stimidx + 30 * fs / 1000];     % samples for clipping waveform from entire trace
    case 'stp'
        % 5 pulses of 500 us at 50 Hz. starts after 10 ms. recording length
        % 200 ms. repeated once every 30 s
%         tstim = 1 / 50 * 1000;      % ms
%         nstim = 5;
        baseline = [1 floor(10 * fs / 1000) - dt];
%         stimidx = [10 : tstim : (nstim - 1) * tstim + 10]; 
%         stimidx  = ceil(stimidx / (1 / fs * 1000));
%         ampidx = [];
%         wvidx = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1 : length(filenames)
    
    % load data
    [~, basename] = fileparts(filenames{i});
    lfp = getLFP('basepath', basepath, 'basename', basename,...
        'extension', 'wcp', 'forceL', true, 'fs', fs, 'saveVar', false,...
        'ch', 1, 'cf', 0, 'concat', true, 'dc', false);   
    lfp.data = rmDC(lfp.data, 'dim', 1, 'win', baseline);
    
    % manually inspect and remove unwanted traces
    if inspect
        [lfp.data, rm{i}] = rmTraces(lfp.data, 'x', lfp.timestamps,...
            'ylim', ylimit);
        lfp.data(:, rm{i}) = [];
    else
        rm{i} = [];
    end
    
    % analyze data; this is your part lior :)
    trace{i} = lfp.data;
%     traceavg(i, :) = mean(lfp.data, 2);
%     wv{i} = lfp.data(wvidx);
%     wvavg(i, :) = mean(lfp.data(wvidx, :), 2);
    
%     switch protocol
%         case 'io'
%             amp(i) = range(traceavg(i, ampidx));
%             ampcell{i} = range(trace{i}(ampidx, :));
%         case 'stp'
%             for ii = 1 : nstim
%                 s1 = stimidx(ii) + dt;
%                 s2 = floor(min([stimidx(ii) + tstim / 1000 * fs - dt, length(trace{i})]));
%                 amp(i, ii) = range(traceavg(i, s1 : s2));
%                 ampcell{i}(ii, :) = range(trace{i}(s1 : s2, :));
%             end
%     end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arrange struct and save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch protocol
    case 'io'
        tstamps = lfp.timestamps*1000-30;
    case 'stp'
        tstamps = lfp.timestamps*1000-10;
end
intensOrig = intens;
[intens, ia] = sort(intens);

% arrange struct - old

% fepsp.origFiles = filenames;
% fepsp.intens = intens;
% fepsp.fs = lfp.fs;
% fepsp.trace = trace(:, ia);
% fepsp.traceavg = traceavg(ia, :);
% fepsp.wv = wv(:, ia);
% fepsp.wvavg = wvavg(ia, :);
% fepsp.tstamps = [0 : 1 / fs : 0.03] * 1000;       
% fepsp.amp = amp(ia);
% fepsp.ampmat = cell2nanmat(ampcell(ia));
% fepsp.rm = rm;
% fepsp.stimidx = stimidx;
% fepsp.ampidx = ampidx;
% fepsp.wvidx = wvidx;

% arrange struct - new (to Match with fEPSPfromDat)
fepsp.info.basename = filenames;
fepsp.info.recSystem = 'wcp';
fepsp.info.protocol = protocol;
fepsp.info.fsOrig = lfp.fs_orig;
fepsp.info.fs = fs;
fepsp.info.spkgrp = {1};
fepsp.info.stimIdx = []; %While we can put the stimIDx here of the traces,
%it won't match the stimIDx in fEPSPfromDat - there it is the index of stimuli in the long Dat file.
fepsp.info.stimTs = [];  % for stp, however same as above - stimTs are the index in the long Dat.
fepsp.info.stamps = [];             % to recover tdt
fepsp.info.intensOrig = intensOrig;
fepsp.info.rm = rm; %Order matching original
fepsp.info.lowPass = 450;
fepsp.info.inspect = inspect;
fepsp.intens = intens;
fepsp.tstamps = tstamps;
fepsp.traces = trace(:, ia);

if saveVar && ~anaflag
    save(fepspname, 'fepsp');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if graphics
%     figure
%     suptitle(expName)
%     subplot(2, 1, 1)    % waveform
%     ph = plot(fepsp.tstamps, fepsp.wvavg, 'k');
%     alphaIdx = [0.3 : 0.5 / length(ph) : 0.8];
%     for i = 1 : length(ph)
%         ph(i).Color(4) = alphaIdx(i);
%     end
%     legend(split(num2str(fepsp.intens)))
%     xlabel('Time [ms]')
%     ylabel('Voltage [mV]')
%     box off
%     ylim(ylimit)
% 
%     subplot(2, 1, 2)    % amplitude
%     switch protocol
%         case 'io'
%             boxplot(fepsp.ampmat, 'PlotStyle', 'traditional')
%             bh = findobj(gca, 'Tag', 'Box');
%             alphaIdx = fliplr(alphaIdx);
%             for i = 1 : length(bh)
%                 patch(get(bh(i), 'XData'), get(bh(i), 'YData'),...
%                     'k', 'FaceAlpha', alphaIdx(i))
%             end
%             xticklabels(split(num2str(fepsp.intens)))
%             xlabel('Intensity [uA]')
%         case 'stp'
%             for i = 1 : length(filenames)
%                 plot(amp(i, :) / amp(i, 1))
%                 hold on
%                 ylim([0 1])
%                 xlabel('Stim. No.')
%             end
%     end
%     ylabel('Amplidute [mV]')
%     box off
%     
%     if saveFig
%         export_fig(expName, '-tif', '-r300', '-transparent')
%     end

if anaflag
    AnalyseParm = {'fepsp', fepsp,'saveFig',saveFig,'saveVar',saveVar,'savename',fepspname, AddAnalyseParm{:}};
    fepsp = fEPSP_analysis(AnalyseParm{:});
end

end


%% Old output description for history keeping
%   fepsp           struct with fields:
%       origFiles   cell with raw data file names
%       intens      vec of intensities (one per file)
%       fs          sampling frequency
%       tstamps     timestamps for clipped response [ms]
%       trace       array of recorded trace (150 ms for io). one cell
%                   per file, each cell is mat of samples x traces
%       traceavg    mat of average trace (files x samples)
%       wv          array of clipped trace (40 ms for io). one cell per
%                   file, each cell is mat of samples x traces
%       wvavg       mat of average clipped trace (files x samples)
%       amp         vec of average amplitudes (one per file)
%       ampmat      mat of amplitudes (trace x file)
%       rm          array of traces removed from recording
%       stimidx     index of stimulus [samples]
%       ampidx      range of indices for amplitude calculation
%       wvidx       range of indices for clipping trace