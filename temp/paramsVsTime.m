function paramsVsTime(data, varargin)

% INPUT:
%   data        numeric mat m x n, columns ploted as separate lines

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addOptional(p, 'yLimit', [0 1], @isnumeric);
addOptional(p, 'lineClr', '', @ischar);
addOptional(p, 'shadeIdx', [], @isnumeric);
addOptional(p, 'avgFlag', false, @islogical);
addOptional(p, 'xLabels', {}, @iscell);

parse(p, varargin{:})
yLimit = p.Results.yLimit;
lineClr = p.Results.lineClr;
shadeIdx = p.Results.shadeIdx;
avgFlag = p.Results.avgFlag;
xLabels = p.Results.xLabels;

hp = plot(data);
if size(data, 2) == length(lineClr)
    for i = 1 : size(data, 2)
        hp(i).Color = lineClr(i);
        hp(i).Color(4) = 1;
    end
end
hold on
if ~isempty(shadeIdx)
    fill([shadeIdx fliplr(shadeIdx)]', sort([yLimit yLimit]), 'y',...
        'FaceAlpha', 0.2, 'EdgeAlpha', 0)
end
if avgFlag
    stdshade(data', 0.1, 'k');
end
ylim(yLimit)
x = xlim;
xlim([1 x(2)])
box off
set(gca, 'TickLength', [0 0])

end