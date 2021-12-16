function assignVars(varArray, isession)

% assigns vars to base workspace.

% assign from varArray with isession as index (see getSessionVars
% for more info). varArray is a cell of variables and isession
% is an index to the relevant session.


if iscell(varArray) && isnumeric(isession)
    
    if isempty(varArray{isession, 1})
        assignin('caller', 'session', [])
    else
        assignin('caller', 'session', varArray{isession, 1}.session)
    end
    if isempty(varArray{isession, 2})
        assignin('caller', 'cm', [])
    else
        assignin('caller', 'cm', varArray{isession, 2}.cell_metrics)
    end
    if isempty(varArray{isession, 3})
        assignin('caller', 'spikes', [])
    else
        assignin('caller', 'spikes', varArray{isession, 3}.spikes)
    end
    if isempty(varArray{isession, 4})
        assignin('caller', 'fr', [])
    else
        assignin('caller', 'fr', varArray{isession, 4}.fr)
    end
    if isempty(varArray{isession, 5})
        assignin('caller', 'datInfo', [])
    else
        assignin('caller', 'datInfo', varArray{isession, 5}.datInfo)
    end
    if isempty(varArray{isession, 6})
        assignin('caller', 'ss', [])
    else
        assignin('caller', 'ss', varArray{isession, 6}.ss)
    end
    if isempty(varArray{isession, 7})
        assignin('caller', 'sr', [])
    else
        assignin('caller', 'sr', varArray{isession, 7}.fr)
    end
    if isempty(varArray{isession, 8})
        assignin('caller', 'st', [])
    else
        assignin('caller', 'st', varArray{isession, 8}.st)
    end
    if isempty(varArray{isession, 9})
        assignin('caller', 'ripp', [])
    else
        assignin('caller', 'ripp', varArray{isession, 9}.ripp)
    end
end