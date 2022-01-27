function setMatlabGraphics(factoryFlag)

if isempty(factoryFlag)
    factoryFlag = true;
end

if ~factoryFlag
    
    set(groot, 'DefaultTextInterpreter', 'none');
    set(groot, 'DefaultAxesTickLabelInterpreter', 'none');
    set(groot, 'DefaultLegendInterpreter', 'none');
    set(groot, 'DefaultAxesFontName', 'Diverda Sans Com')
    set(groot, 'DefaultTextFontName', 'Diverda Sans Com')
    set(groot, 'DefaultAxesFontSize', 13)
    set(groot, 'DefaultAxesOuterPosition', [0 0 1 1])
    set(groot, 'DefaultAxesLabelFontSizeMultiplier', 1.4)
    set(groot, 'DefaultAxesFontWeight', 'normal')
    set(groot, 'DefaultAxesTitleFontWeight', 'bold')
    set(groot, 'DefaultAxesTitleFontSizeMultiplier', 1.6) 
    set(groot, 'DefaultFigurePaperUnits', 'normalize')
    set(groot, 'DefaultFigureUnits', 'normalize')
    set(groot, 'DefaultFigurePaperPosition', [0 0 1 1])
%     set(groot, 'DefaultFigurePosition', [0.2 0.075 0.6 0.85])
    set(groot, 'DefaultFigurePosition', [0.1 0.1 0.8 0.8])
    set(groot, 'DefaultFigurePaperType', 'a4')
    set(groot, 'DefaultFigureColor', [1 1 1])       % figure white instead of gray
    set(groot, 'DefaultFigurePaperOrientation', 'portrait')
    set(groot, 'DefaultLineLineWidth', 0.01)
    
else
    % factory values
    
    set(groot, 'DefaultTextInterpreter', 'remove');
    set(groot, 'DefaultAxesTickLabelInterpreter', 'remove');
    set(groot, 'DefaultLegendInterpreter', 'remove');
    set(groot, 'DefaultAxesFontName', 'remove')
    set(groot, 'DefaultTextFontName', 'remove')
    set(groot, 'DefaultAxesFontSize', 'remove')
    set(groot, 'DefaultAxesOuterPosition', 'remove')
    set(groot, 'DefaultAxesLabelFontSizeMultiplier', 'remove')
    set(groot, 'DefaultAxesFontWeight', 'remove')
    set(groot, 'DefaultAxesTitleFontWeight', 'remove')
    set(groot, 'DefaultAxesTitleFontSizeMultiplier', 'remove')
    set(groot, 'DefaultFigurePaperUnits', 'remove')
    set(groot, 'DefaultFigureUnits', 'remove')
    set(groot, 'DefaultFigurePaperPosition', 'remove')
    set(groot, 'DefaultFigurePosition', 'remove')
    set(groot, 'DefaultFigurePaperType', 'remove')
    set(groot, 'DefaultFigureColor', 'remove')      
    set(groot, 'DefaultFigurePaperOrientation', 'remove')    
    set(groot, 'DefaultLineLineWidth', 'remove')
   
end
end