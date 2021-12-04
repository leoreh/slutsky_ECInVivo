function assignVars(varArray, isession)

% assigns vars to base workspace.

% assign from varArray with isession as index (see getSessionVars
% for more info). varArray is a cell of variables and isession
% is an index to the relevant session.


if iscell(varArray) && isnumeric(isession)
    
    if isempty(varArray{isession, 1})
        assignin('base', 'session', [])
    else
        assignin('base', 'session', varArray{isession, 1}.session)
    end
    if isempty(varArray{isession, 2})
        assignin('base', 'cm', [])
    else
        assignin('base', 'cm', varArray{isession, 2}.cell_metrics)
    end
    if isempty(varArray{isession, 3})
        assignin('base', 'spikes', [])
    else
        assignin('base', 'spikes', varArray{isession, 3}.spikes)
    end
    if isempty(varArray{isession, 4})
        assignin('base', 'fr', [])
    else
        assignin('base', 'fr', varArray{isession, 4}.fr)
    end
    if isempty(varArray{isession, 5})
        assignin('base', 'datInfo', [])
    else
        assignin('base', 'datInfo', varArray{isession, 5}.datInfo)
    end
    if isempty(varArray{isession, 6})
        assignin('base', 'ss', [])
    else
        assignin('base', 'ss', varArray{isession, 6}.ss)
    end
    if isempty(varArray{isession, 7})
        assignin('base', 'sr', [])
    else
        assignin('base', 'sr', varArray{isession, 7}.fr)
    end
    if isempty(varArray{isession, 8})
        assignin('base', 'st', [])
    else
        assignin('base', 'st', varArray{isession, 8}.st)
    end
    if isempty(varArray{isession, 9})
        assignin('base', 'ripp', [])
    else
        assignin('base', 'ripp', varArray{isession, 9}.ripp)
    end
end