function f = addLns(varargin)

% adds lines and labels to current figure
%
% INPUT
%   lns         vector of indices to lines
%   lbs         cell of line labels
%   ax          lines on x-axis {'x'} or y-axis ['y']
%
% 23 dec 19 LH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addOptional(p, 'lns', []);
addOptional(p, 'lbs', {});
addOptional(p, 'ax', 'x', @isstr);

parse(p, varargin{:})
lns = p.Results.lns;
lbs = p.Results.lbs;
ax = p.Results.ax;


if length(lns) ~= length(lbs) 
    error('length of lns and lbs must be equal')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f = gcf;
nlines = length(lns);
l = repmat(lns, 2, 1);

if strcmp('x', ax)
    b = repmat(ylim, nlines, 1);
    plot(l, b' ,'--k', 'HandleVisibility','off')
    if ~isempty(lbs)
        text(l(1, :), b(1, 2)*ones(nlines, 1), lbs)
    end
elseif strcmp('y', ax)
    b = repmat(xlim, nlines, 1);
    plot(b', l ,'--k', 'HandleVisibility','off')
    if ~isempty(lbs)
        text(b(1, 2)*ones(nlines, 1), l(1, :), lbs)
    end
end

end

% EOF