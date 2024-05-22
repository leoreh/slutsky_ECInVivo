function plot_stg(stg, varargin)

% plots the stg and unit data for specified (or random) pairs of significant synapses.
%
% INPUT:
%   stg         struct containing stg data
%   basepath    string specifying the file directory {pwd}
%   nSyn        Number of random synapse pairs to plot, for E and I each
%   synIdx      m x 2 mat specifying pairs of units to plot (overrides nSyn)
%   saveFig     logical flag indicating whether to save the figure to disk
% 
% DEPENDENCIES
%   plot_ccg
% 
% TO DO lIST
%   convert to gui
%
% 10 may 24 LH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addParameter(p, 'basepath', pwd, @ischar);
addParameter(p, 'nSyn', 10, @isnumeric);
addParameter(p, 'synIdx', [], @isnumeric);
addParameter(p, 'saveFig', false, @islogical);

parse(p, varargin{:});
basepath = p.Results.basepath;
nSyn = p.Results.nSyn;
synIdx = p.Results.synIdx;
saveFig = p.Results.saveFig;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% determine synapses to plot
if isempty(synIdx)
     
    % randomly select E and I synapses 
    eIdx = stg.eIdx(randperm(size(stg.eIdx, 1), min(nSyn, size(stg.eIdx, 1))), :);
    iIdx = stg.iIdx(randperm(size(stg.iIdx, 1), min(nSyn, size(stg.iIdx, 1))), :);
        
    % combine selected synapses
    syn2plot = [eIdx; iIdx];

else
    syn2plot = synIdx;  % direct use of provided synapse indices 
end
nSyn = size(syn2plot, 1);

% extract general params from stg struct
roi_binIdx = stg.info.roi_binIdx;
cc50bins = stg.cc.cc50bins;
cc150bins = stg.cc.cc150bins;
cc20bins = stg.cc.cc20bins;

% load other variables
vars = ["swv_metrics"; "units"; "session"];
varsName = ["swv"; "units"; "session"];
[v, basepaths] = getSessionVars('basepaths', {basepath}, 'varsFile', vars,...
    'varsName', varsName);
fs = v.session.extracellular.sr;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plotting each synapse pair
for isyn = 1 : nSyn

    % extract specific synapse data
    u1 = syn2plot(isyn, 1);
    u2 = syn2plot(isyn, 2);   
    dccc = stg.cc.dccc(:, u1, u2);
    pred = stg.info.pred(:, u1, u2);
    eBins = stg.eBins(:, u1, u2);
    iBins = stg.iBins(:, u1, u2);
    
    % determine synapse E or I, where I takes precedence if there are both
    % E and I significant gins
    if stg.iSig(u1, u2)
        synType = -1;
    elseif stg.eSig(u1, u2)
        synType = 1; 
    else
        error('No significant STG between unit %d and %d.', u1, u2);
    end
    
    % determine unit type from cell explorer
    if find(v.units.clean(:, u1)) == 1
        clr(1) = 'b';
    elseif find(v.units.clean(:, u1)) == 2
        clr(1) = 'r';
    else
        clr(1) = 'k';
    end
    if find(v.units.clean(:, u2)) == 1
        clr(2) = 'b';
    elseif find(v.units.clean(:, u2)) == 2
        clr(2) = 'r';
    else
        clr(2) = 'k';
    end

    % Open figure
    fh = figure;
    set(fh, 'WindowState', 'maximized');
    th = tiledlayout(2, 4);
    th.TileSpacing = 'tight';
    th.Padding = 'none';
    title(th, string(sprintf('%s\nSynapse %d-%d',...
        stg.info.basename, u1, u2)), 'interpreter', 'none', 'FontSize', 20);
    set(fh, 'DefaultAxesFontSize', 16);

    % acg unit 1
    axh = nexttile(th, 1); cla; hold on
    plot_ccg(stg.cc.cc150(:, u1, u1), cc150bins, 'clr', clr(1),...
        'pred', [], 'sigbins1', [], 'sigbins2', []);
    title(sprintf('Unit #%d (Presynaptic)', u1))

    % acg unit 2
    axh = nexttile(th, 4); cla; hold on
    plot_ccg(stg.cc.cc150(:, u2, u2), cc150bins, 'clr', clr(2),...
        'pred', [], 'sigbins1', [], 'sigbins2', []);
    title(sprintf('Unit #%d (Postsynaptic)', u2))

    % ccg 50
    axh = nexttile(th, 2); cla; hold on
    plot_ccg(stg.cc.cc50(:, u1, u2), cc50bins, 'clr', 'k',...
        'pred', [], 'sigbins1', roi_binIdx(eBins), 'sigbins2', roi_binIdx(iBins));

    % deconvoluted cc 
    axh = nexttile(th, 3); cla; hold on
    plot_ccg(dccc, cc50bins, 'clr', 'k', 'pred', pred,...
        'sigbins1', roi_binIdx(eBins), 'sigbins2', roi_binIdx(iBins))
    if synType == -1
        title(sprintf('iSTG = %.4f', stg.iStg(u1, u2)))
    elseif synType == 1
        title(sprintf('eSTG = %.4f', stg.eStg(u1, u2)))
    end
    
    % ccg 150
    axh = nexttile(th, 6); cla; hold on
    plot_ccg(stg.cc.cc150(:, u1, u2), cc150bins, 'clr', 'k',...
        'pred', [], 'sigbins1', [], 'sigbins2', []);
    
    % ccg 120
    axh = nexttile(th, 7); cla; hold on
    plot_ccg(stg.cc.cc20(:, u1, u2), cc20bins, 'clr', 'k',...
        'pred', [], 'sigbins1', [], 'sigbins2', []);

    % waveforms
    if ~isempty(v.swv)

        % waveform unit 1
        axh = nexttile(th, 5); cla; hold on
        wv = v.swv.wv(u1, :);
        wv_std = v.swv.wv_std(u1, :);
        x_val = [1 : length(wv)] / fs * 1000;
        plot(x_val, wv, clr(1), 'LineWidth', 2)
        patch([x_val, flip(x_val)], [wv + wv_std, flip(wv - wv_std)],...
            clr(1), 'EdgeColor', 'none', 'FaceAlpha', .2, 'HitTest', 'off')
        xlabel('Time [ms]')
        ylabel('Voltage [mV]')

        % waveform 2
        axh = nexttile(th, 8); cla; hold on
        wv = v.swv.wv(u2, :);
        wv_std = v.swv.wv_std(u2, :);
        plot(x_val, wv, clr(2), 'LineWidth', 2)
        patch([x_val, flip(x_val)], [wv + wv_std, flip(wv - wv_std)],...
            clr(2), 'EdgeColor', 'none', 'FaceAlpha', .2, 'HitTest', 'off')
        xlabel('Time [ms]')
        ylabel('Voltage [mV]')

    end
    
    % save figure
    drawnow
    if saveFig
        figpath = fullfile(basepath, 'graphics', 'stg');
        mkdir(figpath);
        figname = fullfile(figpath, sprintf('stg_%d_%d.fig', u1, u2));
        savefig(fh, figname, 'compact')
    end
end

end
