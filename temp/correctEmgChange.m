
sSig = load([basename, '.sleep_sig.mat']);

timepnt = datInfo.nsamps / fs;

sec1 = sSig.emg_rms(1 : timepnt(1) - 1);
sec2 = sSig.emg_rms(timepnt(1) : end);

prct = 97;
rng = [prctile(sec2, 100 - prct), prctile(sec2, prct)];

sec1 = bz_NormToRange(sec1, rng);


% sSig.emg_rms(timepnt(1) : end) = sec2;

sSig.emg_rms(1 : timepnt(1) - 1) = sec1;

AccuSleep_viewer(sSig, [], [])



save([basename, '.sleep_sig.mat'], '-struct', "sSig")