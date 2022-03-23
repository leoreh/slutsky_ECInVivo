

% fepsp_sessions

% experiment folder
exppath = 'G:\Data\hb_a2';
cd(exppath)
basepaths = dir('*_io*');
basepaths = fullfile({basepaths.folder}, {basepaths.name});
nfiles = length(basepaths);

% stim session names
[~, basenames] = fileparts(basepaths);
ids = cellfun(@(x) x(end - 2 : end), basenames, 'uni', false);

% load data
varsFile = ["fepsp_traces"; "fepsp_results"; "lfp"];
varsName = ["traces"; "results"; "lfp"];
vIO = getSessionVars('basepaths', basepaths, 'varsFile', varsFile,...
    'varsName', varsName);

% load stp
cd(exppath)
basepaths = dir('*_stp*');
basepaths = fullfile({basepaths.folder}, {basepaths.name});
nfiles = length(basepaths);
varsFile = ["fepsp_traces"; "fepsp_results"; "lfp"];
varsName = ["traces"; "results"; "lfp"];
vSTP = getSessionVars('basepaths', basepaths, 'varsFile', varsFile,...
    'varsName', varsName);

% graphics
fh = figure;
th = tiledlayout(3, 2);

% raw io traces
nexttile
xval = [1 : size(vIO(end - 2).traces{2}, 1)] / 10000 * 1000;
plot(xval, mean(vIO(end - 2).traces{2}, 2, 'omitnan'))
hold on
plot(xval, mean(vIO(end).traces{2}, 2, 'omitnan'))
xlabel('Time [ms]')
ylabel('Amplitude [mV]')
legend({'Baseline', 'Baclofen'})
title('Traces @ 70 uA')

% io curve
nexttile
xval = [50, 70, 90];
plot(xval, [vIO(end - 2).results.avg_traces.Amp{:}])
hold on
plot(xval, [vIO(end).results.avg_traces.Amp{:}])
xticks([30, 50, 70, 90])
xlabel('Intensity [uA]')
ylabel('Amplitude [mV]')
title('IO Curve')

% ------- stp
% load

% stp traces
nexttile
xval = [1 : size(vSTP(end - 2).traces{1}, 1)] / 10000 * 1000;
plot(xval, mean(vSTP(end - 2).traces{1}, 2, 'omitnan'))
hold on
plot(xval, mean(vSTP(end).traces{1}, 2, 'omitnan'))
axis tight
xlabel('Time [ms]')
ylabel('Amplitude [mV]')
legend({'Baseline', 'Baclofen'})
title('Traces @ 50 uA')

% stp facilitation curve
nexttile
xval = [50, 70, 90];
yval = [vSTP(end - 2).results.avg_traces.Amp{:}];
plot(xval, yval ./ yval(1) * 100)
hold on
yval = [vSTP(end).results.avg_traces.Amp{1}];
plot(xval, yval / yval(1) * 100)
xticks([30, 50, 70, 90])
xlabel('Stim No.')
ylabel('Norm. Amplitude [%]')
title('Facilitation')

% free run
cd(exppath)
basepaths = dir('*_freerun*');
basepaths = fullfile({basepaths.folder}, {basepaths.name});
basepath = basepaths{1};
cd(basepath)
load("220316_freerun.lfp.mat")

nexttile(th, [1 2])
[s, tstamps, freq] = specBand('basepath', basepath, 'sig', lfp.data,...
    'fs', lfp.fs, 'graphics', true, 'logfreq', false);
yLimit = ylim;
hold on
tidx = [cumsum(lfp.filelength) / 60 / 60]';
tidx = [1; tidx(1 : end - 1)];
plot([tidx, tidx], ylim, '--k', 'LineWidth', 2)
hold on
text(tidx(end - 2), yLimit(2) + 0.01 * yLimit(2), "Baclofen")









exppath = 'G:\Data\hb_a2';
basename = 'test_stp';
fepsp_protocol = 'stp';
wcpfiles = [24];
intens = [50];

% load single file
fepsp_wcpPipeline('basepath', exppath, 'wcpfiles', wcpfiles,...
    'intens', intens, 'recname', basename,...
    'fepsp_protocol', fepsp_protocol, 'saveFlag', true)





% 
% % experiment folder
% exppath = 'G:\Data\hb_a2';
% cd(exppath)
% basepaths = dir('*_stp*');
% basepaths = fullfile({basepaths.folder}, {basepaths.name});
% nfiles = length(basepaths);
% 
% for ifile = 1 : nfiles
%     
%     ifile = 6;
%     cd(basepaths{ifile})
%     [~, basename] = fileparts(basepaths{ifile});
%     load([basename, '.lfp.mat'])
%     [~, wcpfiles] = fileparts(lfp.files)
%     intens = [50, 70];
% 
%     if ~iscell(wcpfiles)
%         wcpfiles = {wcpfiles};
%     end
% 
%     fepsp_wcpPipeline('basepath', exppath, 'wcpfiles', wcpfiles,...
%         'intens', intens, 'fsOut', [], 'recname', basename,...
%         'fepsp_protocol', lfp.fepsp_protocol)
% 
% end
