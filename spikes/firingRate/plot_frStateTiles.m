function [mdl, y2] = plot_frStateTiles(varargin)

% state-dependent regulation of mfr. to be continued

% INPUT
%   basepath        recording session {pwd}
%   stateMfr        2 x n numeric of mfr per unit in two states (rows).
%                   row 1 should be aw. if empty will load from fr
%                   struct.
%   fr              struct. see calc_fr. if empty will load from basepath
%   units           2 x n logical mat of indices to stateMfr. row 1 is RS
%                   and row 2 if FS. if empty, all units in stateMfr will
%                   be plotted.
%   ntiles          numeric. number of tiles to divide the units {6}.
%   stateNames      2 x 1 cell of char. if empty will assume wake and nrem.
%   saveFig         logical {true}
%   graphics        logical {true}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addOptional(p, 'basepath', pwd);
addOptional(p, 'ntiles', 4, @isnumeric);
addOptional(p, 'saveFig', true, @islogical);
addOptional(p, 'graphics', true, @islogical);

parse(p, varargin{:})
basepath    = p.Results.basepath;
ntiles      = p.Results.ntiles;
saveFig     = p.Results.saveFig;
graphics    = p.Results.graphics;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load vars
varsFile = ["fr"; "units"];
varsName = ["fr"; "units"];
[v, ~] = getSessionVars('basepaths', {basepath}, 'varsFile', varsFile,...
    'varsName', varsName);
[~, basename] = fileparts(basepath);

% select params
idxUnits = v.units.rs;
unitType = 'RS';
unitClr = 'k';
stateNames = {'WAKE', 'NREM'};

% organize fr in states and ratio
x = v.fr.states.mfr(idxUnits, 1);  % mfr in aw
y = v.fr.states.mfr(idxUnits, 4);  % mfr in nrem
stateRat = squeeze(v.fr.states.ratio(4, 1, :));
stateMfr = [x y];

% divide units into percentiles and oragnize state ratio accordingly
[~, tileIdx, tiles] = vec2tileMat(v.fr.mfr, ntiles, v.units.rs);
frTiles = cellfun(@(x) stateRat(x), tileIdx, 'uni', false);
dataMat = cell2nanmat(frTiles, 2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% linear regression (on log-log mfr)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% equality line
eqLine = 10 .^ [floor(log10(min(stateMfr(stateMfr(:) ~= 0)))),...
    ceil(log10(max(stateMfr(:))))];

% ignore zeros 
idxInf = isinf(log10(x)) | isinf(log10(y));

mdl = fitlm(log10(x(~idxInf)), log10(y(~idxInf)));
fitCoeff = flipud(mdl.Coefficients.Estimate);
x2 = logspace(-4, 4, 9);
y2 = 10 .^ [polyval(fitCoeff, log10(x2))];
fitLine = [x2; y2];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if graphics
    setMatlabGraphics(true)

    fh = figure;
    subplot(2, 1, 1)
    plot_boxMean('dataMat', dataMat, 'clr', unitClr, 'allPnts', true)
    hold on
    plot(xlim, [0, 0], '--k')
    ylabel({sprintf('%s - %s /', stateNames{2}, stateNames{1}),...
        sprintf('%s + %s', stateNames{2}, stateNames{1})})
    xlabel('Pecentiles')
    ylim([-100, 100])
    title(sprintf('%s, %s = %d units', basename, unitType, sum(idxUnits)))

    subplot(2, 1, 2)
    plot(x, y, '.', 'Color', unitClr, 'MarkerSize', 10)
    hold on
    plot(eqLine, eqLine, 'k')
    plot(fitLine(1, :), fitLine(2, :), '--', 'Color', unitClr)
    plot([tiles(2 : end - 1); tiles(2 : end - 1)], eqLine, '--g')
    set(gca, 'YScale', 'log', 'XScale', 'log')
    xlabel('WAKE firing rate [Hz]')
    ylabel('NREM firing rate [Hz]')
    xlim(eqLine)
    ylim(eqLine)
    title(sprintf('Slope = %.2f; R2 = %.2f',...
        mdl.Coefficients.Estimate(2), mdl.Rsquared.Ordinary))

    if saveFig
        figpath = fullfile('graphics');
        mkdir(figpath)
        figname = fullfile(figpath, sprintf('%s_frStateTiles', basename));
        export_fig(figname, '-tif', '-transparent', '-r300')
    end
end

end

% EOF