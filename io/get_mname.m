function mname = get_mname(basepaths)

% returns the mouse names from an array of basepaths

mname = cell(size(basepaths));

for ipath = 1 : length(basepaths)
    
    % Extract mouse name using regular expression
    % Pattern explanation: Match 'lh' followed by one or more digits
    matches = regexp(basepaths{ipath}, 'lh\d+', 'match');
    
    % Store the first match found (assuming there's only one per string)
    if ~isempty(matches)
        mname{ipath} = matches{1};
    else
        % Fallback: extract the last folder name from the path
        pathParts = strsplit(basepaths{ipath}, '\');
        mname{ipath} = pathParts{end};
    end
end

end
