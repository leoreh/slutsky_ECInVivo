function [inb, b] = selectCluster(data, varargin)

% allows user to select points by drawing borders around a cluster in plot
% 
% INPUT
%   data        2D matrix to plot
%   sep         2 element vector used to plot linear separatix. for example:
%               sep = [-2 2.5] will plot y = -2 * data(:, 1) + 2.5
%   newplot     logical. plot data or use current figure
% 
% OUTPUT
%   inb         logical vector where 1 = in border and 0 = outside border
%   b           border coordinates
% 
% DEPENDENCIES
%   PointInput      
%
% 08 apr 19 LH. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addOptional(p, 'sep', [], @isnumeric);
addOptional(p, 'h', [], @ishandle);

parse(p, varargin{:})
sep = p.Results.sep;
h = p.Results.h;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(h)
    figure;
    hold on
    plot(data(1, :), data(2, :), 'k.');
end

if ~isempty(sep)
    xlim([0, max(data(1, :)) + 0.1])
    ylim([0, max(data(2, :)) + 0.1])
    xb = get(gca, 'XLim');
    yb = get(gca, 'YLim');
    plot(xb, [sep(1) * xb(1) + sep(2), sep(1) * xb(2) + sep(2)])
end

linepos = 1;
b = [];    % matrix of polyline coords
h_lastline = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% user selection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
zoomoff = 1;
while 1
    [x, y, button]= PointInput(1);
    res = button * zoomoff;
    title({'Select cluster',...
        'left click to draw boundary',...
        'right click to complete',...
        'middle click to end'});    
    switch res
        case 1 % left button
            b(linepos, :) = [x y];
            if linepos > 1
                h_lastline(end+1) = line([b(linepos-1,1),b(linepos,1)],[b(linepos-1,2),b(linepos,2)]);
                set(h_lastline(end), 'Color', 'k');
            end
            linepos = linepos + 1;
        case 2  % middle button
            break
        case 3  % right button
            if linepos > 2
                b(linepos, :) = b(1, :);
                h_lastline(end + 1) = line([b(linepos-1,1),b(linepos,1)],[b(linepos-1,2),b(linepos,2)]);
                set(h_lastline(end),'Color','k');
            end            
    end
end
inb = inpolygon(data(1, :), data(2, :), b(:, 1),b(:, 2));

end

