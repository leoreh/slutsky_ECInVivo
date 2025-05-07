

% manually inspect .lfp files

% create eeg and spec signals that are same length and params as lfp

% manually inspect spec across experiment, per mouse, eeg and lfp top and
% bottom

% sleep scoring

% create separate psd struct (for emg only, probably)


%%% mouse comments
% lh99 - unclear why forfieted mouse. 
% lh105 - only lfp, eeg, and emg. mouse recorded in dtd.
% lh106 - emg and eeg in separate binary file. mouse also includes lfp.
% recorded in tdt

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate psd per file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% grab mice with eeg signals
grps = [mcu_sessions('eeg')];

% state params
ftarget = [1 : 0.5 : 100];

% iterate
for igrp = 1 : length(grps)

    queryStr = grps{igrp};
    basepaths = mcu_sessions(queryStr);
    nfiles = length(basepaths);
    mpath = fileparts(basepaths{1});

    for ifile = 1 : nfiles
        
        % session params
        session = CE_sessionTemplate(pwd, 'viaGUI', false,...
            'forceDef', false, 'forceL', false, 'saveVar', true);      
        basepath = session.general.basePath;
        nchans = session.extracellular.nChannels;
        [~, basename] = fileparts(basepath);

        % files
        cd(basepath)
        sigfile = fullfile(basepath, [basename, '.sleep_sig.mat']);

        % load lfp signal
        sig = load(sigfile, 'eeg');
        sig = sig.eeg;
        load(sigfile, 'fs');



        % calc psd accordign to emg state separation
        psdEmg = psd_states('basepath', basepath, 'sstates', [1, 2],...
            'sig', sig, 'fs', fs, 'saveVar', true,...
            'graphics', true, 'forceA', true, 'ftarget', ftarget,...
            'flgEmg', true);
    end
end
