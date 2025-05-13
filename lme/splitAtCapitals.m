function firstPart = splitAtCapitals(inputStr)
% SPLITATCAPITALS Splits a string at capital letters and returns the first part
%   Example: splitAtCapitals('BurstRate') returns 'Burst'
%
% INPUT:
%   inputStr - A string to split
%
% OUTPUT:
%   firstPart - The first part of the string before the first capital letter

% Find the index of the first capital letter after the first character
capIdx = regexp(inputStr(2:end), '[A-Z]');
if isempty(capIdx)
    firstPart = inputStr; % If no capital letters found, return the whole string
else
    firstPart = inputStr(1:capIdx(1)); % Return everything up to the first capital
end

end 