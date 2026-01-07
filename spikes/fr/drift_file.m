function drft = drift_file(varargin)

% wrapper for drift_calc. repeats the calculation for RS and FS units, and
% for sleep states. relies on fr and units files in basepath
%
% INPUT:
%   basepath        string. path to recording folder {pwd}
%   boutTimes       cell of n x 2 mats. each cell describes the bouts of
%                   a state (s). if empty will try to load from
%                   sleep_states.mat and will only analyze aw and nrem
%   graphics        logical. plot {false}
%
% DEPENDENCIES:
%   drift_calc
%   drift_plot
%
% TO DO LIST:
%
% 22 may 24 LH      based on Lee's code. see also geva (ziv) et al., Neuron, 2013

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addOptional(p, 'basepath', pwd);
addOptional(p, 'boutTimes', []);
addOptional(p, 'graphics', false, @islogical);

parse(p, varargin{:})
basepath        = p.Results.basepath;
boutTimes     = p.Results.boutTimes;
graphics        = p.Results.graphics;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~, basename] = fileparts(basepath);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% analysis params
winsize = 3600;
thrLin = 0;
thrFr = 0.005;
thrWin = 0;
timeAlt = 1;

% load data
% load data
varsFile = ["fr"; "units"; "sleep_states"];
v = basepaths2vars('basepaths', {basepath}, 'vars', varsFile);
if isfield(v, 'SleepState'), [v.ss] = v.SleepState; v = rmfield(v, 'SleepState'); end

% state params
sstates = [1, 4];       % selected states for calculating drift
cfg = as_loadConfig;
snames = cfg.names(sstates);

if isempty(boutTimes)
    if ~isempty(v.ss)
        boutTimes = v.ss.boutTimes(sstates);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% analyze
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iunit = 1 : 2
    unitIdx = v.units.clean(iunit, :);

    % calc drift across recording
    drft(iunit, 1) = drift_calc(v.fr.strd(unitIdx, :), v.fr.tstamps,...
        'graphics', false, 'winsize', winsize, 'thrLin', thrLin,...
        'thrFr', thrFr, 'thrWin', thrWin);

    % calc drift per state
    for istate = 1 : length(sstates)

        stateIdx = InIntervals(v.fr.tstamps, boutTimes{istate});
        fr_mat = v.fr.strd(unitIdx, stateIdx);

        switch timeAlt
            case 1
                % ALT 1: concatenate bouts and ignore actual time when
                % calculating drift
                tstamps = [30 : 60 : size(fr_mat, 2) * 60];

            case 2
                % ALT 2: maintain original tstamps such that drift includes
                % elpased time
                tstamps = v.fr.tstamps(stateIdx);

            case 3
                % ALT 3: grab fr in states from fr struct. same as ALT 2
                % but more accurate. overrides the boutTimes input
                fr_mat = v.fr.states.fr{sstates(istate)}(unitIdx, :);
                tstamps = v.fr.states.tstamps{sstates(istate)};
        end

        drft(iunit, istate + 1) = drift_calc(fr_mat, tstamps, 'graphics', false,...
            'winsize', winsize, 'thrLin', thrLin, 'thrFr', thrFr, 'thrWin', thrWin);
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot drift per unit and state

ulabel = ["RS"; "FS"];
d = catfields(drft, 'addim');
yLimit = [min(d.m_corr, [], 'all'), 1];

if graphics

    setMatlabGraphics(true)
    fh = figure;
    set(fh, 'WindowState', 'maximized');
    tlayout = [1, length(sstates) + 1];
    th = tiledlayout(tlayout(1), tlayout(2));
    th.TileSpacing = 'tight';
    th.Padding = 'none';
    set(fh, 'DefaultAxesFontSize', 16);
    title(th, basename, 'interpreter', 'none', 'FontSize', 20);

    for iunit = 1
        tbias = (iunit - 1) * tlayout(2) + 1;
        axh = nexttile(th, tbias, [1, 1]); cla; hold on
        drift_plot(drft(iunit, 1), axh);
        ylabel(sprintf('%s Drift Rate [1 / h]', ulabel(iunit)))
        title(axh, 'Full Recording')
        ylim(yLimit)

        for istate = 1 : length(sstates)
            tbias = (iunit - 1) * tlayout(2) + 1 + istate;
            axh = nexttile(th, tbias, [1, 1]); cla; hold on
            drift_plot(drft(iunit, istate + 1), axh);
            ylabel(sprintf('%s Drift Rate [1 / h]', ulabel(iunit)))
            title(axh, snames(istate))
            ylim(yLimit)

        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% organize data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for each numeric field, that two last dimesions are unit and then state
for istate = 1 : length(sstates) + 1
    s(istate) = catfields(drft(:, istate), 'addim');
end
s = catfields(s, 'addim');
flds = fieldnames(s);
for ifld = 1 : length(flds)
    fld = flds{ifld};
    if isnumeric(s.(fld)) || islogical(s.(fld))
        s.(fld) = squeeze(s.(fld));
    end
end
drft = s;

end

% EOF

