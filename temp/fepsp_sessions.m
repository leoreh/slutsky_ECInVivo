


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% single session analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

basepath = 'G:\Data\lh103\220329';

wcpfiles = [16, 17];
intens = [38, 36];

% io
fepsp_wcpPipeline('basepath', basepath, 'wcpfiles', wcpfiles,...
    'fepsp_protocol', 'stp', 'recname', '220329_stp1', 'intens', intens,...
    'saveFlag', true)

% freerun
fepsp_wcpPipeline('basepath', basepath, 'wcpfiles', wcpfiles,...
    'fepsp_protocol', 'freerun', 'recname', 'freerun', 'intens', [],...
    'saveFlag', true)



% plot spec ---------------------------------------------------------------
% load
basepath = pwd;
[~, basename] = fileparts(basepath);
specfile = fullfile(basepath, [basename, '.spec.mat']);
load(specfile)

% graphics
setMatlabGraphics(false)
fh = figure;

plot_spec(spec, false, false)
hold on
axis tight
yLimit = ylim;
tidx = [cumsum(lfp.filelength) / 60 / 60]';
tidx = [0; tidx(1 : end - 1)];
plot([tidx, tidx], ylim, '--k', 'LineWidth', 1)
hold on
% text(tidx, ones(1, length(tidx)) * (yLimit(2) + 0.05 * yLimit(2)), string(wcpfiles))

% save
figpath = fullfile(basepath, 'graphics');
mkdir(figpath)
figname = fullfile(figpath, sprintf('%s_spec', basename));
export_fig(figname, '-tif', '-transparent', '-r300')

% psd timebins ------------------------------------------------------------

fs = lfp.fs;
nsec = cumsum(lfp.filelength);
clear winCalc
winCalc{1} = [1, nsec(1)];
winCalc{2} = [nsec(1), nsec(3)];
winCalc{3} = [nsec(3), nsec(4)];
% winCalc{4} = [nsec(6), nsec(7)];
winCalc{4} = [nsec(8), nsec(9)];
winCalc{5} = [nsec(9), nsec(10)];

winCalc = cellfun(@round, winCalc, 'uni', false);

[psdBins, faxis] = psd_timebins('sig', lfp.data, 'fs', fs,...
    'winCalc', winCalc, 'graphics', true);
legend

psdBins ./ sum(psdBins, 2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% multi-session analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load data ---------------------------------------------------------------

% experiment folder
exppath = 'G:\Data\lh103\220329';
cd(exppath)
basepaths = dir('*_stp*');
basepaths = fullfile({basepaths.folder}, {basepaths.name});
nfiles = length(basepaths);

% stim session names
[~, basenames] = fileparts(basepaths);
ids = cellfun(@(x) x(end - 2 : end), basenames, 'uni', false);

% load data
varsFile = ["fepsp_traces"; "fepsp_results"; "lfp"];
varsName = ["traces"; "results"; "lfp"];
v = getSessionVars('basepaths', basepaths, 'varsFile', varsFile,...
    'varsName', varsName);


% to prism ----------------------------------------------------------------

for ifile = 1 : nfiles
    
    % intens
    [intens, sidx] = sort(v(ifile).lfp.intens);

    % mean trace per intensity
    traces = cellfun(@(x) mean(x, 2, 'omitnan'), v(ifile).traces, 'uni', false);
    traces = cell2mat(traces);
    traces = traces(:, sidx);

    % timestamps
    protocol_info = fepsp_getProtocol('protocol_id', v(ifile).lfp.fepsp_protocol);
    tstamps = v(ifile).lfp.timestamps * 1000 - protocol_info.stim_times(1); 
    tstamps = tstamps(1 : length(traces));


    % amp
    amp = cell2nanmat(v(ifile).results.all_traces.Amp, 2)';
    amp = amp(sidx, :);

    % slope
    slope = cell2nanmat(v(ifile).results.all_traces.Slope, 2)';
    slope = slope(sidx, :);
    


    v(ifile).results.avg_traces.Amp{1}
    v(ifile).results.avg_traces.Amp{1} / v(ifile).results.avg_traces.Amp{1}(1)

    v(ifile).results.avg_traces.Slope{1}
    v(ifile).results.avg_traces.Slope{1} / v(ifile).results.avg_traces.Slope{1}(1)


end














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
