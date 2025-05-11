




% sleep scoring

% create separate psd struct (for emg only, probably)


%%% mouse comments
% LH99 - EEG RUINED ON BAC ON BUY MAYBE I CAN STILL USE FOR ACUTE EFFECT
% lh100 - first 12 hr of baseline ruined (eeg only), also ruined after 4
% hours of washout (does not return)
% lh105 - only lfp, eeg, and emg. mouse recorded in tdt.
% lh106 - emg and eeg in separate binary file. mouse also includes lfp.
% recorded in tdt

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create EEG signal and spectrogram
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% grab mice with eeg signals
mNames = [mcu_sessions('eeg')];
nMice = length(mNames);

% import toolbox for filtering
import iosr.dsp.*

% iterate
for iMouse = 1 : nMice

    queryStr = mNames{iMouse};
    basepaths = mcu_sessions(queryStr);
    nfiles = length(basepaths);
    mpath = fileparts(basepaths{1});

    for ifile = 1 : nfiles

        % filenames
        basepath = basepaths{ifile};
        [~, basename] = fileparts(basepath);
        cd(basepath)
        sessionfile = [basename, '.session.mat'];
        ssigfile = [basename, '.sleep_sig.mat'];
        eegfile = [basename, '.sleep_eeg.mat'];
        specfile = [basename, '.spec_eeg.mat'];

        % load
        load(sessionfile)
        load(ssigfile, 'emg')

        % preparations
        fs = session.extracellular.srLfp;
        nchans = session.extracellular.nChannels;

        % load channels and average
        if strcmp(mNames{iMouse}, 'lh106')
            eegName = [basename, '.emg.dat'];
            eegCf = [450];
            eeg = double(bz_LoadBinary(eegName, 'duration', Inf,...
                'frequency', 1250, 'nchannels', 2, 'start', 0,...
                'channels', 2, 'downsample', 1));
            
            

            tstamps_sig = [1 : length(emg)] / fs;

            % interpolate
            eeg = [interp1([1 : length(eeg)] / fs, eeg, tstamps_sig,...
                'spline')]';

            filtRatio = eegCf / (fs / 2);
            eeg = iosr.dsp.sincFilter(eeg, filtRatio);
        else
        
            eegName = [basename, '.lfp'];
            
            ch = [2 : 8];
            dur = Inf;
            start = 0;
            eeg = double(binary_load(eegName, 'duration', dur,...
                'fs', fs, 'nCh', nchans, 'start', start,...
                'ch', ch, 'downsample', 1, 'bit2uv', 0.195));
        end

        if length(eeg) ~= length(emg)
            error('length')
        end

        % calculate spectrogram
        spec_eeg = calc_spec('sig', eeg, 'fs', 1250, 'graphics', true,...
            'saveVar', false, 'padfft', -1, 'winstep', 5,...
            'ftarget', [], 'ch', {1},...
            'force', true, 'logfreq', true);

        % save eeg signal
        save(eegfile, "eeg")
        save(specfile, "spec_eeg")

    end
end

% manually inspect spec across experiment, per mouse, eeg and lfp top and
% bottom
for iMouse = 1 : nMice
    mName = mNames{iMouse};
    [~, info] = sessions_catVarTime('mname', mName,...
        'dataPreset', {'spec', 'spec_eeg', 'hypnogram_emg'}, 'graphics', true, 'dataAlt', 1,...
        'basepaths', {}, 'xTicksBinsize', 6, 'markRecTrans', true,...
        'saveFig', false);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create EMG epochs 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is also the opportunity to get rid of bad episodes in the rocording.
% This seems only necassary for lh100. This means lh100 EEG bouts will no
% longer correspond to LFP bouts. 

% However, due to the removal of outlier bouts, which is done based on the
% PSDs, the population of EMG and LFP bouts are expected to be different
% for each mouse.
% Unless, we apriori use the bouts from the psdEmg file, rather then from
% sleep_statesEmg
% But again, maybe the removal of outliers should not be trusted, given it
% is blind to bout duration

% for lh100, perhaps better to limit existing bouts manually rather than
% recalculating them. 

% just a reference - no use as of now
interDur = 3;
minDur = 5;
for ifile = 1 : nfiles

    % files
    basepath = basepaths{ifile};
    [~, basename] = fileparts(basepath);
    cd(basepath)
    ssEmg = as_emg('basepath', basepath, 'flgInspct', false,...
        'interDur', interDur, 'minDur', minDur);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate EEG PSD per bout 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% grab mice with eeg signals
mNames = [mcu_sessions('eeg')];
nMice = length(mNames);

% params
ftarget = [1 : 0.5 : 100];

% iterate
for iMouse = 4 : nMice

    queryStr = mNames{iMouse};
    basepaths = mcu_sessions(queryStr);
    nfiles = length(basepaths);
    mpath = fileparts(basepaths{1});

    for ifile = 1 : nfiles

        % files
        basepath = basepaths{ifile};
        [~, basename] = fileparts(basepath);
        cd(basepath)
        sigfile = fullfile(basepath, [basename, '.sleep_eeg.mat']);
        statesfile = [basename '.sleep_statesEmg.mat'];
        psdfile = [basename '.psdEMG.mat'];
        psdEEGfile = [basename '.psdEMG_eeg.mat'];
        % load lfp signal
        eeg = load(sigfile, 'eeg');
        fs = 1250;
        
        % get bout times from emg states
        load(statesfile)
        if strcmp(queryStr, 'lh100')
            btimes = [];
        else
            % btimes = psd.bouts.times;
            btimes = ssEmg.bouts.times;
        end

        % calc psd accordign to emg state separation
        psd = psd_states('basepath', basepath, 'sstates', [1, 2],...
            'sig', eeg.eeg, 'fs', fs, 'saveVar', false,...
            'graphics', false, 'forceA', true, 'ftarget', ftarget,...
            'flgEmg', true, 'btimes', btimes, 'flgOtl', true);
        
        % save
        save(psdEEGfile, "psd");

    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FOOOF Analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


flgEeg = false;
mNames = [mcu_sessions('eeg')];
mNames = [mcu_sessions('wt'), mcu_sessions('mcu')];

if flgEeg
    varPsd = '.psdEMG_eeg.mat';
    varF1f = '.f1f_eeg.mat';
else
    varPsd = '.psdEMG.mat';
    varF1f = '.f1f.mat';
end

fRange = [1 100];
nMice = length(mNames);

for iMouse = 1 : nMice
    mName = mNames{iMouse};
    mPaths = mcu_sessions(mName);
    nFiles = length(mPaths);

    for ifile = 1 : nFiles
        basepath = mPaths{ifile};
        [~, basename] = fileparts(basepath);
        original_path = pwd;
        cd(basepath);
        
        psdFile = fullfile(basepath, [basename, varPsd]);
        f1ffile = fullfile(basepath, [basename, varF1f]);

        % extract psd struct
        psd = load(psdFile);
        flds = fieldnames(psd);
        psd = psd.(flds{1});

        freqs = psd.info.faxis;
        nStates = size(psd.psd, 1);
        
        clear f1f f1fBout f1fState
        for iState = 1 : nStates
            psdState = psd.bouts.psd{iState};
            stateBtimes = psd.bouts.times{iState};

            [boutLbls, boutStats] = bouts_separate(stateBtimes,...
                'sepMet', 'mcshane', 'flgGraphics', false);
            boutMasks = {~boutLbls, boutLbls};
            nBoutGrps = length(boutMasks);
            % nBoutGrps = 1;

            for iBoutGrp = 1 : nBoutGrps 
                psdAvg = mean(psdState(boutMasks{iBoutGrp}, :), 1);
                
                f1fBout(iBoutGrp) = fooof_calc(psdAvg, freqs, ...
                                  'f_range', fRange, ...
                                  'flg_plot', false, ...
                                  'saveVar', false); 
            end
            f1fState(iState) = catfields(f1fBout, 1);
        end
        f1f = catfields(f1fState, 'addim', [], [3, 1, 2], true);
           
        save(f1ffile, 'f1f', '-v7.3');
        cd(original_path); 
    end 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FOOOF Inspect
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

flgEeg = false;
mGrp = 'wt';            % 'wt' / 'mcu' / 'eeg'
mNames = mcu_sessions(mGrp);
nMice = length(mNames);

if flgEeg
    vars = ["f1f_eeg"];
else
    vars = ["f1f"];
end

% assumes origianl struct organized as        [state x boutGrp x band]
% after day concatenation will be, per mouse  [day x state x boutGrp x band]
% after mouse concatenation will be           [mouse x day x state x boutGrp x band]

clear f1f f1fMouse
idx_rmDays = [2, 6];              % remove bac on and off
for iMouse = 1 : nMice
    mName = mNames{iMouse};
    mPaths = mcu_sessions(mName);
    mPaths(idx_rmDays) = [];
    nFiles = length(mPaths);
    
    v = basepaths2vars('basepaths', mPaths, 'vars', vars);    
    f1fMouse(iMouse) = catfields([v(:).f1f], 'addim', [], [4, 1, 2, 3]);

end
f1f = catfields(f1fMouse, 'addim', [], [5, 1, 2, 3, 4]);

% [dimAvg, dimGrp, dimRow]
size(f1f.psd_orig)
fh = fooof_plotGrp(f1f, 1, 2, 3);


% [dimAvg, dimGrp, dimRow]
size(f1f.psd_orig)
fh = fooof_plotGrp(f1f, 1, 4, 3);

% don't forget lh100. maybe this means bouts needs to be saved within f1f
% struct (separate eeg from lfp)
% find a way to track boutStats per file. 
% perhaps also solve spec outliers.
% they do seem meaningful and bleed to nearby timebins. maybe recalculate
% in 2-3s windows. 

% one struct for handeling states. start only with emg. 
% original separation.
% spec outliers - create new bouts
% in psd_states, update this states structure


% grab only during specific bouts
size(f1f.psd_orig)
psdWt = squeeze(f1f.psd_orig(:, 1, 1, 1, :));
freqs = squeeze(f1f.freqs(1, 1, 1, 1, :));
figure
plot(freqs, psdWt);
axh = gca;
set(axh, 'XScale', 'log', 'YScale', 'log')
legend

psdWt = squeeze(f1f.psd_orig(:, 1, 2, :, :));
psdWtData = squeeze(psdWt(:, 1, :));
freqs = squeeze(f1f.freqs(1, 1, 1, 1, :));

fh = figure;
axh = subplot(1, 1, 1);
ph = plot_stdShade('dataMat', psdWtData, 'xVal', freqs, 'alpha', 0.3,...
    'axh', axh);
set(axh, 'XScale', 'log', 'YScale', 'log')
hold on
ph = plot_stdShade('dataMat', psdMcu, 'xVal', freqs, 'alpha', 0.3,...
    'axh', axh);


% grab only during specific bouts
size(f1f.psd_orig)
psdMcu = squeeze(f1f.psd_orig(:, 1, 2, 1, :));

psdWt = squeeze(f1f.psd_orig(:, 1, 2, :, :));
psdWtData = squeeze(psdWt(:, 1, :));
freqs = squeeze(f1f.freqs(1, 1, 1, 1, :));

fh = figure;
axh = subplot(1, 1, 1);
ph = plot_stdShade('dataMat', psdWtData, 'xVal', freqs, 'alpha', 0.3,...
    'axh', axh);
set(axh, 'XScale', 'log', 'YScale', 'log')
hold on
ph = plot_stdShade('dataMat', psdMcu, 'xVal', freqs, 'alpha', 0.3,...
    'axh', axh);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LME on FOOOF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% organzie all bands in lme_tbl

% load data for each group
grps = {'wt', 'mcu'};             % 'wt' / 'mcu' / 'eeg'
grps = {'eeg'};                   % 'wt' / 'mcu' / 'eeg'
idx_rmDays = [2, 6];              % remove bac on and off
clear grppaths
for iGrp = 1 : length(grps)

    mNames = mcu_sessions(grps(iGrp)); 
    for iMouse = 1 : length(mNames)
        basepaths = mcu_sessions(mNames{iMouse});
        basepaths(idx_rmDays) = [];
        grppaths{iGrp}(iMouse, :) = string(basepaths)';
    end
end

% FOOOF ~ Group * Day + (1|Mouse)
% organize for lme1
frml = 'FOOOF ~ State * Day + (1|Mouse)';
[lme_tbl, lme_cfg] = lme_org('grppaths', grppaths, 'frml', frml,...
    'flg_emg', true, 'var_field', 'bands', 'vCell', {});

% select state
% istate = "High EMG";
% plot_tbl = lme_tbl(lme_tbl.State == categorical(istate), :);

% normalize
plot_tbl = lme_normalize('lme_tbl', lme_tbl, 'normVar', 'Day',...
    'groupVars', {'Group', 'State'});
% plot_tbl = lme_tbl;

% run lme
contrasts = 'all';
[lme_results, lme_cfg] = lme_analyse(plot_tbl, lme_cfg, 'contrasts', contrasts);

% plot
fh = lme_plot(plot_tbl, lme_cfg.mdl, 'ptype', 'bar');
