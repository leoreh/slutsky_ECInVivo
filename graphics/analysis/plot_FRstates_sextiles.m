function plot_FRstates_sextiles(varargin)

% receives a mat of mfr per unit in different states and plots the state
% ratio according to sextiles and the correlation between mfr in states.
% sextiles are created according to the mfr in state 1 (1st row of
% stateMfr), typically wake. 

% INPUT
%   basepath        recording session {pwd}
%   stateMfr        2 x n numeric of mfr per unit in two states (rows). 
%                   row 1 should be wake.
%   units           2 x n logical mat of indices to stateMfr. row 1 is RS 
%                   and row 2 if FS. if empty, all units in stateMfr will
%                   be plotted. 
%   ntiles          numeric. number of tiles to divide the units {6}.
%   stateNames      2 x 1 cell of char. if empty will assume wake and nrem.
%   saveFig         logical {true}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addOptional(p, 'basepath', pwd);
addOptional(p, 'stateMfr', [], @isnumeric);
addOptional(p, 'units', []);
addOptional(p, 'ntiles', 6, @isnumeric);
addOptional(p, 'stateNames', {}, @iscell);
addOptional(p, 'saveFig', true, @islogical);

parse(p, varargin{:})
basepath    = p.Results.basepath;
stateMfr    = p.Results.stateMfr;
units       = p.Results.units;
ntiles      = p.Results.ntiles;
stateNames  = p.Results.stateNames;
saveFig     = p.Results.saveFig;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate state ratio
stateRat = [(stateMfr(2, :) - stateMfr(1, :)) ./...
    (stateMfr(2, :) + stateMfr(1, :))]';

% organize units
if isempty(units)
    units = logical(ones(1, length(stateMfr)));
    unitsAny = units;
end

if size(units, 1) == 2
    unitType = {'RS', 'FS'};
    unitClr = {'b', 'r'};
    unitsAny = (units(1, :) | units(2, :));
else
    unitType = {'Units'};
    unitClr = {'k'};
end

if isempty(stateNames)
  stateNames = {'WAKE', 'NREM'}; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
setMatlabGraphics(false)

eqLine = 10 .^ [floor(log10(min(stateMfr(stateMfr(:) ~= 0)))),...
    ceil(log10(max(stateMfr(:))))];
yLimit = [-max(abs(stateRat(unitsAny))), max(abs(stateRat(unitsAny)))];

fh = figure;
for iunit = 1 : size(units, 1)
    
    subplot(2, size(units, 1), iunit)
    
    x = stateMfr(1, units(iunit, :));
    y = stateMfr(2, units(iunit, :));
    
    % sextiles. note calculation of tiles includes all units if they were
    % not excluded during selectUnits.
    tiles = prctile(x, [1 / ntiles : 1 / ntiles : 1 - 1 / ntiles] * 100);
    tiles = [0, tiles, Inf];

    clear dataMat
    for itile = 1 : ntiles
        unitsTile = units(iunit, :) & stateMfr(1, :) >= tiles(itile) &...
            stateMfr(1  , :) < tiles(itile + 1);
        dataMat{itile} = stateRat(unitsTile);
    end
    dataMat = cell2nanmat(dataMat);
    plot_boxMean('dataMat', dataMat, 'clr', unitClr{iunit}, 'allPnts', true)
    hold on
    plot(xlim, [0, 0], '--k')
    ylabel({sprintf('%s - %s /', stateNames{2}, stateNames{1}),...
        sprintf('%s + %s', stateNames{2}, stateNames{1})})
    xlabel(sprintf('Sixtiles (by %s)', stateNames{1}))
    ylim(yLimit)
    subtitle(sprintf('%s = %d', unitType{iunit}, sum(units(iunit, :))))
    
    subplot(2, size(units, 1), iunit + size(units, 1))
    infidx = isinf(log10(x));   
    mdl = fitlm(log10(x(~infidx)), log10(y(~infidx)));
    fitCoeff = flipud(mdl.Coefficients.Estimate); 
    y2 = 10 .^ [polyval(fitCoeff, log10(eqLine(1) : eqLine(2)))];
    plot(x, y, '.', 'Color', unitClr{iunit}, 'MarkerSize', 10)
    hold on
    plot(eqLine, eqLine, 'k')
    plot([eqLine(1) : eqLine(2)], y2, '--', 'Color', unitClr{iunit})
    plot([tiles(2 : end - 1); tiles(2 : end - 1)], eqLine, '--g')
    set(gca, 'YScale', 'log', 'XScale', 'log')
    xlabel('WAKE firing rate [Hz]')
    ylabel('NREM firing rate [Hz]')
    xlim(eqLine)
    ylim(eqLine)
    subtitle(sprintf('Slope = %.2f; R2 = %.2f',...
        mdl.Coefficients.Estimate(2), mdl.Rsquared.Ordinary))
end

end

% EOF