
% prep
basepath = 'D:\Data\fepsp_test\lh85_211012_113400';
[~, basename] = fileparts(basepath);
cd(basepath)
load('lh85_catted.mat')

protocol_id = 'io';
intens = IntensIO;
data_in = CattedfepspIO;
stimIdx = StimIndexIO;

protocol_id = 'stp';
intens = IntensSTP;
data_in = CattedfepspSTP;
stimIdx = StimIndexSTP;




% figure
stimN = 2;
% plot(CattedfepspIO)
% hold on
% yLimit = ylim;
% plot([StimIndexIO{stimIdx}; StimIndexIO{stimIdx}], yLimit, '--k')

% step 1
traces = fepsp_org2traces('data_in', data_in, 'fs', fs,...
    'protocol_id', protocol_id, 'stim_locs', stimIdx);

figure, plot(traces{stimN})

% step 2
marking_win = fepsp_markings("traces", traces, "fs", fs,...
    "protocol_id", protocol_id,"base_path", basepath,...
    "intens", intens, "traces_Xlimit", [], "traces_Ylimit", [],...
    "dt", 2, "max_time_tol", 3, "fast_mark", true);

load([basename, '_fepsp_markings.mat'])

% step 3
Results = fepsp_analyse("traces", traces, "fs", fs,...
    "protocol_id", protocol_id, "marked_peaks", marked_peaks,...
    "marked_starts", marked_starts, "Avg_marked_peaks", Avg_marked_peaks,...
    "Avg_marked_starts", Avg_marked_starts,...
    "base_path", basepath, "save_var", true, "slope_area", [0.2 0.9]);

% step 4
analysed_fepsp = fepsp_show("traces", traces, "fs", fs,...
    "protocol_id", protocol_id, "Avg_marked_peaks", Avg_marked_peaks,...
    "Avg_marked_starts", Avg_marked_starts, "Results", Results,...
    "base_path", basepath, "intens", intens);

