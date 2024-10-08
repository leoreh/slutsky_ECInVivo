% IED_wrapper

basepath = pwd;
cd(basepath)
[~, basename] = fileparts(basepath);

% load lfp from WCP/lfp
channels_in_file = 4;
lfp_chan2load = 2;
eeg_chan2load = 3;
emg_chan2load = 4; % note this is optional, leave empty ( [] ) & comment relevant lines to skip emg


file_type = 'lfp';
% lfp = getLFP('basepath', basepath, 'nchans', channels_in_file, 'ch', [lfp_chan2load, emg_chan2load], 'chavg', {},...
%     'fs', 1250, 'interval', [0 inf], 'extension', file_type,...
%     'savevar', false,'concat',true, 'forceL', false, 'basename', basename,'pli', false); 
% fs = lfp.fs;
% sig = lfp.data(:,1);
% emg = lfp.data(:,2);

% If you use .lfp file, its fs is already 1250, and you don't care for the
% 450 low pass / have already a filtered file:
% You can use the following lines to get better performance & saving
% smaller file size.
% Comment "getLFP" and the lines next to it, and uncomment the following lines:
file_name = fullfile(basepath, join([basename ".lfp"],""));
map_format = map_creator(file_name,channels_in_file,"int16");
maped_file = memmapfile(file_name,"Format",[map_format 'Mapped']);
sig = memmap_row(maped_file,lfp_chan2load);
emg = memmap_row(maped_file,emg_chan2load); % note assumption that emg is same file. if it isn't, just do the full process with a new file.
eeg = memmap_row(maped_file,eeg_chan2load);
fs = 1250;

thr_mv = 0;
thr_Z  = 5;
ied = IED.detect(sig, fs, "LFP",...
    eeg, fs, "EEG", ...
    emg, fs, "EMG", ...
    "sig2use", "LFP", "thr",[thr_Z, thr_mv],"thrDir","both");
app = reduct_displayer(ied,...
    "reduction_z_norm",false,...
    "reduction_function","t-SNE","tsne_Algorithm",'barneshut',"tsne_dist_fun",'correlation',...
    "cluster_method","dbscan"...
    );
clusters2use = [];
ied.accepted = ismember(app.cluster_idx, clusters2use);
ied = IED.curate_only_init_accepted(ied,"saveVar",true,"basepath",basepath,"basename",basename);

% % The ied, from class IED.data, is the heart of everything - use it to
% % stop & restart working:
% % To restart from last point, simply run the line with IED.curate again
% % (make sure you have everything you need in workspace).
% % IED.curate with the above syntax will save your "ied" that have
% % everything you need automatically. However you can & should save it manually
% % if needed, by:
% ied.file_loc = fullfile(basepath,join([basename "ied.mat"],"."));
% save(ied.file_loc,"ied")
%
% % to load it, simply load the mat file:
% mat2load = fullfile(basepath,join([basename "ied.mat"],"."));
% load(mat2load,"ied")
%
% % Note that you can see what stage of analysis you are, by looking at: 
% ied.status

% analyze after curation:
bin_size = ied.fs*60; % bin for rate, in [samples]
smooth_win = 7;       % size of movemean window, in [bins]
clip_marg = 0.05;     % how much time to show around the event
ied = IED.analyze(ied,"binsize",bin_size,"smf",smooth_win,"marg",clip_marg,...
    "saveFig",true,"saveVar",true);
% Use right click to enlarge any of the axes.
% Use right click on a IED to change it to declined - note that after doing
% this, you will need to re-run IED.analyze.

% EOF