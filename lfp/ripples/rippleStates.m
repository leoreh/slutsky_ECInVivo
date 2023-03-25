function ripp = rippleStates(ripp, varargin)

% analysis ripples in relation to vigilance states. relies on the struct ss
% from ss_classify.m. plots the results
%
% INPUT:
%   basepath            path to recording {pwd}
%   graphics            logical. plot graphics {true} or not (false)
%   saveVar             logical. save variables (update ripp)
%
% OUTPUT:
%   ripp            struct
%
% DEPENDENCIES:
%   none
%
% 12 jan 23 LH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addOptional(p, 'basepath', pwd);
addOptional(p, 'graphics', true, @islogical);
addOptional(p, 'saveVar', true, @islogical);

parse(p, varargin{:})
basepath            = p.Results.basepath;
graphics            = p.Results.graphics;
saveVar             = p.Results.saveVar;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% files
cd(basepath)
[~, basename] = fileparts(basepath);
rippfile = fullfile(basepath, [basename, '.ripp.mat']);

% load
if isempty(ripp)
    if exist(rippfile, 'file')
        load(rippfile, 'ripp')
    else
        error('ripp file missing')
    end
end

% load session vars
varsFile = ["session"; "sleep_states"; "spikes.cellinfo"; "spktimes"];
varsName = ["session"; "ss"; "spikes"; "spktimes"];
v = getSessionVars('basepaths', {basepath}, 'varsFile', varsFile,...
    'varsName', varsName);
ss = v.ss;

% state params
cfg = as_loadConfig([]);
sstates = [1 : 5];          % selected states

% params
fsSpk       = v.session.extracellular.sr;
fs          = ripp.info.fs;
epochs      = ripp.epochs;
nepochs     = size(epochs, 1);
recWin      = ripp.info.recWin;
nbinsMap    = size(ripp.maps.freq, 2);
durWin      = [-75 75] / 1000;

% initialize
ripp.states = [];
ripp.states.stateNames = ss.info.names;
nstates = length(ss.stateEpochs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ripple during states
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for istate = sstates

    % limit state epochs to recWin
    epochIdx = InIntervals(ss.stateEpochs{istate}, recWin);
    if ~isempty(ss.stateEpochs{istate})

        % rate in states
        [ripp.states.rate{istate}, ripp.states.binedges{istate},...
            ripp.states.tstamps{istate}] =...
            times2rate(ripp.peakPos, 'binsize', ripp.info.binsizeRate,...
            'winCalc', ss.stateEpochs{istate}(epochIdx, :), 'c2r', true);

        % idx of rippels in state
        ripp.states.idx{istate} =...
            InIntervals(ripp.peakPos, ss.stateEpochs{istate}(epochIdx, :));
    else
        ripp.states.idx{istate} = nan;
        ripp.states.rate{istate} = nan;
        ripp.states.binedges{istate} = nan;
        ripp.states.tstamps{istate} = nan;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save and graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if saveVar
    save(rippfile, 'ripp')
end

if graphics

    setMatlabGraphics(true)
    fh = figure;

    % rate
    sb1 = subplot(1, 4, [1, 2]);
    plot(ripp.rate.tstamps / 60 / 60, ripp.rate.rate, 'k')
    xlabel('Time [h]')
    ylabel('Ripple Rate [Hz]')
    hold on
    for istate = sstates
        ph = plot(ripp.states.tstamps{istate} / 60 / 60,...
            ripp.states.rate{istate}, '.', 'MarkerSize', 10);
        ph.Color = cfg.colors{istate};
    end
    xlabel('Time [h]')
    ylabel('Ripple Rate [Hz]')

    % percent rippels in state
    sb2 = subplot(1, 4, [3, 4]);
    prct_states = sum(cell2nanmat(ripp.states.idx, 2), 1, 'omitnan');
    pie(prct_states, ones(1, length(prct_states)))
    hold on
    ph = findobj(sb2, 'Type', 'Patch');
    set(ph, {'FaceColor'}, flipud(cfg.colors(sstates)))
    legend({ripp.states.stateNames{sstates}}, 'FontSize', 10,...
        'Units', 'normalized', 'Location', 'best',...
        'NumColumns', 2);

    sgtitle(basename)

    % save figure
    figpath = fullfile(basepath, 'graphics');
    mkdir(figpath)
    figname = fullfile(figpath, sprintf('%s_rippleStates', basename));
    export_fig(figname, '-tif', '-transparent', '-r300')

end

end

% EOF