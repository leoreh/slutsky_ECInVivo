function assignVars(varArray, isession, vars)

% assign to workspace vars from varArray with isession as index (see
% getSessionVars for more info). 

if nargin < 3 || isempty(vars)
    vars = ["session";...
        "cm";...
        "spikes";...
        "fr";...
        "datInfo";...
        "ss";...
        "sr";...
        "st";...
        "ripp"];
end

for ivar = 1 : length(vars)
    % fix spkrate so can be w/ firing in same workspace
    if strcmp(vars{ivar}, 'sr')
        assignin('caller', vars{ivar}, varArray{isession, ivar}.fr)
    elseif strcmp(vars{ivar}, 'cm')
        if isempty(varArray{isession, ivar})
            assignin('caller', vars{ivar}, [])
        else
            assignin('caller', vars{ivar}, varArray{isession, ivar}.cell_metrics)
        end
    else
        if isempty(varArray{isession, ivar})
            assignin('caller', vars{ivar}, [])
        else
            assignin('caller', vars{ivar}, varArray{isession, ivar}.(vars{ivar}))
        end
    end
end

