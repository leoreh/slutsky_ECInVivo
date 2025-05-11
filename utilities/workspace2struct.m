function sOut = workspace2struct(sIn)
% WORKSPACE2STRUCT Populates struct fields from caller's workspace variables.
%
% SUMMARY:
%   Populates fields of an input struct `sIn` with values from variables
%   of the same name in the caller's workspace. If a variable is not
%   found in the caller's workspace, a warning is issued, and the
%   original value in `sIn` is retained.
%
% INPUT:
%   sIn         (Required) Input struct.
%
% OUTPUT:
%   sOut        Output struct with fields populated from caller's workspace.
%
% lh 251105

sOut = sIn;
fieldNames = fieldnames(sIn);

% Get the name of the input struct variable for more informative warnings.
structNameInCaller = inputname(1);
if isempty(structNameInCaller)
    structNameInCaller = 'inputStruct'; % Default if called with an unnamed expression
end

for iField = 1:length(fieldNames)
    fieldName = fieldNames{iField};

    % Check if variable exists in the caller's workspace
    if evalin('caller', ['exist(''' fieldName ''', ''var'')'])
        % Assign variable's value from caller's workspace to output struct
        sOut.(fieldName) = evalin('caller', fieldName);
    else
        % Warn if variable not found; field retains original value
        warning('WORKSPACE2STRUCT:notFound', ...
            'Variable ''%s'' not found in caller''s workspace for struct ''%s''. Field ''%s'' retains its original value.', ...
            fieldName, structNameInCaller, fieldName);
    end
end

end     % EOF