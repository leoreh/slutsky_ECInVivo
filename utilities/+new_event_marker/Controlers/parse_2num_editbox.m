function in_parsed = parse_2num_editbox(user_in)
% parse incoming text to extract what numbers it include,
% allowing simple eqations (i.e. that include '+','-', '^', '*'
% and/or '/' ) & scientific notation using 'e'.
% expect incoming numbers to be seperated by whitespaces.

% check how many elements to expect
nElm = numel(split(user_in,whitespacePattern));

% collect input and transfer 2 num, include equations * & ^ & /
user_in = regexp(user_in,'*?([+-]?\d+[\.*^/e\d+]?|[+-]?inf)*','match');

% if numebr do not match, return nan
if numel(user_in) ~= nElm
    in_parsed = nan;
    return
end

% convert to numbers, return them sorted
in_parsed = cellfun(@(x) str2num(x,"Evaluation","restricted"),user_in,'un',0); %#ok<ST2NM> % I want to eval equations
in_parsed = sort([in_parsed{:}]);
end