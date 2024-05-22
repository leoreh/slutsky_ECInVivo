
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% state (emg) epochs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for each session (and mouse), calculate EMG states
% calculate fr in timebins

% gen parans
saveVar = false;
graphics = true;

% state params
minDur = 10;
interDur = 3;
ftarget = [0.5 : 0.5 : 100];

% files
mname = 'lh134';
basepaths = [mcu_sessions(mname)];
nfiles = length(basepaths);

for ifile = 1 : nfiles

    basepath = basepaths{ifile};
    [~, basename] = fileparts(basepath);
    cd(basepath)


    [stateEpochs, epochStats] = as_epochs('minDur', minDur, 'interDur',...
        interDur, 'flgEmg', true, 'graphics', true, 'nbins', 2);
    
    % calc psd according to emg state separation
    psd = psd_states('basepath', basepath, 'sstates', [1, 2],...
        'sig', [], 'fs', fs, 'saveVar', saveVar,...
        'graphics', true, 'forceA', true, 'ftarget', ftarget,...
        'emgThr', [], 'flgEmg', true, 'stateEpochs', stateEpochs);

    fh = figure;
    plot_hypnogram('stateEpochs', psd.epochs.stateEpochs, 'clr', psd.info.clr, 'axh', gca,...
    'sstates', [1 : 2])

end
