function exlTbl = lme2exl(lme)
    % Extract fixed effects and their statistics from a linear mixed-effects model (lme)
    
    % Get coefficient names and statistics
    coeffs = fixedEffects(lme);
    stats = lme.Coefficients;
    
    % Prepare the data for export
    coeffNames = stats.Name;
    estimates = stats.Estimate;
    SE = stats.SE;
    tStat = stats.tStat;
    pValues = stats.pValue;

    
    % Create a table for easy export
    exlTbl = table(coeffNames, estimates, SE, tStat, pValues, ...
        'VariableNames', {'Coefficient', 'Estimate', 'SE', 'tStat', 'pValue'});
    
    % % Write to clipboard (for quick paste into Excel)
    % writetable(exlTbl, 'clipboard', 'Delimiter', '\t');
    % 
    % % Display the table in MATLAB
    % disp(exlTbl);
    fprintf('Results copied to clipboard. You can now paste directly into Excel.\n');
end