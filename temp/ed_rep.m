
% Load
basepath = pwd;
basename = '220824_0906';
basename = '220615_0801';
basename = '220611_0750';

sSig = load([basename, '.sleep_sig.mat']);
% load([basename, '.sleep_labelsMan.mat'])

% Inspect states
AccuSleep_viewer(sSig, [], [])

% ED detection
thr_Z = 7;
ied = IED.detect_move_z(sSig.eeg, sSig.fs, "LFP",...
    sSig.eeg, sSig.fs, "EEG", ...
    sSig.emg, sSig.fs, "EMG", ...
    "sig2use", "LFP", "thr", [thr_Z], "thrDir", "both");

% drop detections coinciding with high EMG 
ied = IED.reject_emg(ied, "thrZ", 3);

app = reduct_displayer(ied, {},...
    "reduction_z_norm", false,...
    "reduction_function", "t-SNE", "tsne_Algorithm", 'barneshut', "tsne_dist_fun", 'correlation',...
    "cluster_method", "dbscan"...
    );
clusters2use = [1];
ied.accepted = ismember(app.cluster_idx, clusters2use);

ied = IED.curate_only_init_accepted(ied, "saveVar", true, "basepath", basepath, "basename", basename);
