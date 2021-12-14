function f = plotCCG(varargin)

% plot cross correlogram. if vector than plots single ccg, if matrix [t x
% ngroups x ngroups] than plots all pairs
%
% INPUT
%   ccg         vector (single ccg) or matrix (pairs). see CCG.m
%   t           vector of time lags [s].
%   c           color of hist {k}. if scalar all will be that color. if
%               vector, must be equal length to nunits and ACC (diagonal) will be that
%               color.
%   u           vector of unit identification (e.g. spikes.UID)
%   saveFig     logical. save figure {1} or not (0)
%   basepath    string. recording session path {pwd}
%
% 23 nov 19 LH.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addOptional(p, 'ccg', []);
addOptional(p, 't', []);
addOptional(p, 'c', {'k'}, @iscell);
addOptional(p, 'u', []);
addOptional(p, 'saveFig', false, @islogical);
addOptional(p, 'basepath', pwd);

parse(p, varargin{:})
ccg = p.Results.ccg;
t = p.Results.t;
c = p.Results.c;
u = p.Results.u;
saveFig = p.Results.saveFig;
basepath = p.Results.basepath;

nunits = size(ccg, 2);
if length(c) == 1
    c = repmat(c, nunits, 1);
end
if isempty(u)
    u = 1 : nunits;
end

t = t * 1000;       % s to ms

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% split units to different figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% x = [1 : nunits];
% for i = 2 : nunits - m + 1
%     y = [i - 1 : nunits];
%     yy = [i : nunits];
%     [q, ~, idx] = intersect(y(1 : m), yy);
%     yy(idx) = [];
%     x = [x, i, yy];
% end
% 
% 
% 1   2   3   4   5
% 1   6   7   8   9  
% 1   10  2   6   7
% 2   8   9   10  3
% 3   6   7   8   9
% 
% 
% n = mod(length(x), m);
% z = reshape(x(1 : end - n), m, floor(length(x) / m)) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isvector(ccg)

    % single ccg
    f = gcf;
    b = bar(t, ccg, 'BarWidth', 1);
    ylabel('Counts')
    xlabel('Time [ms]')
    b.FaceColor = c{1};
    b.EdgeColor = 'none';
    box off
    axis tight
    set(gca, 'TickLength', [0 0])

    % grid of ccg
elseif ndims(ccg) == 3   
    f = figure;
    z = 0;
    for i = nunits : -1 : 1       % row
        for j = 1 : (nunits - z)  % column
            k = j + (z * nunits);
            subplot(nunits, nunits, k)
            b = bar(t, ccg(:, j, i), 'BarWidth', 1);
            if i == j
                b.FaceColor = c{i};
            else
                b.FaceColor = 'k';
            end
            b.EdgeColor = 'none';
            box off
            axis tight
            set(gca, 'TickLength', [0 0])
            if j == 1
                ylabel(num2str(u(i)))
            end
            if j == nunits - z
                xlabel(num2str(u(j)))
            end
            i;
            j;
        end
        z = z + 1;
        z = min([z, nunits - 1]);
    end
end

if saveFig
    filename = 'CCG';
    savePdf(filename, basepath, f)
end

end

% EOF