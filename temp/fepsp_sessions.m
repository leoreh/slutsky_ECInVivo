
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% wrapper for tdt 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% single session analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% wcp
basepath        = 'F:\Data\lh110\220727';
wcpfiles        = [11, 12, 21, 22, 35, 46];
intens          = [62];
recname         = '220727_freerun';
fepsp_wcpPipeline('basepath', basepath, 'wcpfiles', wcpfiles,...
    'fepsp_protocol', 'freerun', 'recname', recname, 'intens', intens,...
    'fsOut', [])

% tdt
basepath        = 'F:\Data\lh110\lh110_220728_093900';
mapch           = [];
rmvch           = [1, 3];
intens          = [300, 500, 700];
blocks          = [1 : 5];
protocol_id     = 'pair';
store           = 'Raw2';
recsuffix       = 'io';
ch              = 2;
fepsp_tdtPipeline('basepath', basepath, 'blocks', blocks,...
    'protocol_id', protocol_id, 'recsuffix', recsuffix, 'intens', intens,...
    'ch', ch, 'mapch', mapch', 'rmvch', rmvch, 'store', store)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% multi-session analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load data ---------------------------------------------------------------

% experiment folder
exppath = 'G:\Data\lh104\lh104_220426_090900';
cd(exppath)
basepaths = dir('*_io*');
basepaths = fullfile({basepaths.folder}, {basepaths.name});

nfiles = length(basepaths);

% stim session names
[~, basenames] = fileparts(basepaths);
ids = cellfun(@(x) x(end - 2 : end), basenames, 'uni', false);

% load data
varsFile = ["fepsp_traces"; "fepsp_results"];
varsName = ["traces"; "results"];
v = getSessionVars('basepaths', basepaths, 'varsFile', varsFile,...
    'varsName', varsName);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% to prism
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ifile = 1 : nfiles
    
    % intens
    [intens, sidx] = sort(v(ifile).results.info.intens);

    % mean trace per intensity
    traces = cellfun(@(x) mean(x, 2, 'omitnan'), v(ifile).traces, 'uni', false);
    traces = cell2mat(traces);
    traces = traces(:, sidx);

    % timestamps
    protocol_info = fepsp_getProtocol('protocol_id', v(ifile).lfp.protocol_id);
    tstamps = v(ifile).lfp.timestamps * 1000 - protocol_info.stim_times(1); 
    tstamps = tstamps(1 : length(traces));

    
    % stp intensity idx
    idx = 1;

    % amp
    amp = cell2nanmat(v(ifile).results.all_traces.Amp(idx), 2)';
    amp = amp(sidx, :);

    % slope
    slope = cell2nanmat(v(ifile).results.all_traces.Slope(idx), 2)';
    slope = slope(sidx, :);
    

    ampNorm = v(ifile).results.avg_traces.Amp{idx} / v(ifile).results.avg_traces.Amp{idx}(1);
    slopeNorm = v(ifile).results.avg_traces.Slope{idx} / v(ifile).results.avg_traces.Slope{idx}(1);


end


% single session
avgtraces = cellfun(@(x) mean(x, 2, 'omitnan'), traces, 'uni', false);
cell2nanmat(avgtraces, 2);

results.avg_traces.Amp{1}
results.avg_traces.Amp{1} / results.avg_traces.Amp{1}(1)

lfp.timestamps * 1000 - 10

cell2nanmat(results.all_traces.Amp, 2)










% params
fileIdx = [1 : 7];
intens = 70;
fileNames = {'Baseline', 'aCSF', 'aCSF2', 'Ket' ,'Ket2', 'ket3', 'ket24'};

% graphics
fh = figure;
th = tiledlayout(3, 2);

% raw io traces
nexttile
xval = [1 : size(v(end - 2).traces{2}, 1)] / 10000 * 1000;
hold on
for ifile = 1 : length(fileIdx)
    [~, intensIdx] = min(abs(intens - (v(fileIdx(ifile)).lfp.intens))); 
    plot(xval, mean(v(fileIdx(ifile)).traces{intensIdx}, 2, 'omitnan'))
end
legend
xlabel('Time [ms]')
ylabel('Amplitude [mV]')
legend(fileNames)
title('Traces @ 70 uA')

% io curve amp
nexttile
hold on
for ifile = 1 : length(fileIdx)
[intensVal, sidx] = sort(v(fileIdx(ifile)).lfp.intens);
    plot(intensVal,...
    [v(fileIdx(ifile)).results.avg_traces.Amp{sidx}])
end
legend(fileNames)
xlabel('Intensity [uA]')
ylabel('Amplitude [mV]')
title('IO Curve')

% io curve slope
nexttile
hold on
for ifile = 1 : length(fileIdx)
[intensVal, sidx] = sort(v(fileIdx(ifile)).lfp.intens);
    plot(intensVal,...
    [v(fileIdx(ifile)).results.avg_traces.Slope{sidx}])
end
legend(fileNames)
xlabel('Intensity [uA]')
ylabel('Slope')
title('IO Curve')




% ------- stp
% load

% load stp
cd(exppath)
basepaths = dir('*_stp*');
basepaths = fullfile({basepaths.folder}, {basepaths.name});
nfiles = length(basepaths);
varsFile = ["fepsp_traces"; "fepsp_results"; "lfp"];
varsName = ["traces"; "results"; "lfp"];
vSTP = getSessionVars('basepaths', basepaths, 'varsFile', varsFile,...
    'varsName', varsName);

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




