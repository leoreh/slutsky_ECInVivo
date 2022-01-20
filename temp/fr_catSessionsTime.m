% fr_catSessionsTime

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mname = 'lh96';
forceL = true;
forceA = true;

pcond = ["tempflag"];
ncond = [""];

% load vars from each session
varsFile = ["fr"; "sr"; "cell_metrics"; "datInfo"; "session"];
varsName = ["fr"; "sr"; "cm"; "datInfo"; "session"];
if ~exist('v', 'var') || forceL
    [v, basepaths] = getSessionVars('mname', mname, 'varsFile', varsFile,...
        'varsName', varsName, 'pcond', pcond, 'ncond', ncond);
end
nsessions = length(basepaths);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% params
fs = v(1).session.extracellular.sr;

% initialize
fr = [];
tstamps = [];
xidx = [];
lastT = 0;
for isession = 1 : nsessions
    fr = [fr, v(isession).sr.strd];
    tstamps = [tstamps, v(isession).sr.tstamps / 60 / 60 + lastT];
    xidx = [xidx, v(isession).datInfo.nsamps / fs / 60 / 60 + lastT];
    lastT = tstamps(end);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
setMatlabGraphics(false)

fh = figure;
plot(tstamps, fr)
hold on
plot([xidx; xidx], ylim, '--k')
legend(split(num2str(1 : 4)))

