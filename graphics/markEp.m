function ep = markEp(xdata, ydata, varargin)

% plots a figure and allows the user to mark episodes. left click marks the
% start of an episode, right click marks the end. middle click terminates
% the function. note that marking occurs with repititions (i.e. marking the
% same region twice is possible)

% INPUT
%   xdata       vector. data for x axis
%   ydata       vector. data for y axis
%   xlab        string {'Time [m]'}. xlabel
%   ylab        string {{Voltage [mV]'}. ylabel
%
% OUTPUT
%   ep          n x 2 mat where n is the number of episodes        
%
% TO DO LIST
%       # aesthetics
%       # handle deletions better
%       # xlim size on a log scale
%       # add subplot of entire recording (zoom marked)
%
% 25 mar 20 LH  updates
% 06 may 20 LH      sort(ep)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addParameter(p, 'xlab', 'Time [m]', @isstring)
addParameter(p, 'ylab', 'Voltage [mV]', @isstring)

parse(p, varargin{:})
xlab = p.Results.xlab;
ylab = p.Results.ylab;

% params 
l = length(xdata);

% initialize vars
ep = [];
prevButton = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fh = figure;
suptitle({'Select episodes',...
    'left = select',...
    'right = complete',...
    'middle = end'});
sb1 = subplot(3, 1, 1 : 2);
plot(xdata, ydata, 'k', 'LineWidth', 1)
hold on
set(gca, 'TickLength', [0 0], 'Color', 'none')
box off
ylabel(ylab)
xlabel(xlab)
axis tight
Y = ylim;

% bottom = center
c = uicontrol('Style', 'slider', 'Position', [20 20 60 20]);
c.Value = 0.5;
c.Callback = @winCenter;
% top = size
b = uicontrol('Style', 'slider', 'Position', [20 60 60 20]);
b.Value = 0.2;
b.Callback = @winCenter;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get user input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

while 1
    [x, ~, button] = ginput(1);
    if button == 1
        if prevButton == 1
            ch = get(gca, 'children');
            delete(ch(1))
        end
        line([x ,x], Y);
        temp(1) = x;
    elseif button == 3
        temp(2) = x;
        if prevButton == 3
            ch = get(gca, 'children');
            delete(ch(1))
        else
            ep = [ep; temp];
        end
        fill([temp fliplr(temp)], [Y(1) Y(1) Y(2) Y(2)]',...
            'k', 'FaceAlpha', 0.3,  'EdgeAlpha', 0);
    elseif button == 2
        break
    end
    prevButton = button;
end
close(fh)

ep = sort(ep);

    % callback function
    function winCenter(src, effect)
        cen = xdata(round(c.Value * l));
        siz = xdata(round(b.Value * l));
        xlim([cen - siz cen + siz])
    end
end


% EOF