% as_changes
% list of changes done to accusleep functions

% general
% removed standardizeSR from all functions
% combined emg and eeg to single struct. 
% pass spectrogram and emg rms instead of raw data

% processEMG
% removed the 20-50 passband filter

% AccuSleep_viewer
% # show frequencies up to 15 (line 115)
% # set(G.A6a,'XTick',[]) => set(G.A6a,'XTick',[]) - line 379
% # set(G.A6, 'YTick', []) => set(G.A6) - line 389
% set graphics to default => setMatlabGraphics(true) - line 27

% AccuSleep_classify
% # output netScores (line 90)
% # remove last state from config file so the net won't be trained on bin

% AccuSleep_train
% # remove last state from config file so the net won't be trained on bin
% # added net info to output (line 247)

% createCalibrationData
% # remove last state from config file so the net won't be trained on bin
% # added tolerance to difference in length of labels and t of spectrogram