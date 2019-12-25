function sep = sepBimodel(varargin)

%  finds the separation point of a bimodel distribution according to the
%  trough in a (log10) histogram. tests if the distribution is lognormal
%  via the lillietest.
% 
% INPUT
%   x           vector to be modled  
%   lognorm     logical. log {1} or linear [0] distribution
%   nbins       scalar. number of bins for histcounts
%   graphics    logical. plot figure {1}.
% 
% OUTPUT
%   sep         trough between distribution
% 
% 22 nov 19 LH.
% 
% TO DO LIST
%   add confidance test for separation
%   find test for bimodality

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addOptional(p, 'x', []);
addOptional(p, 'lognorm', true, @islogical);
addOptional(p, 'nbins', 200, @isnumeric);
addOptional(p, 'graphics', true, @islogical);

parse(p, varargin{:})
x = p.Results.x;
lognorm = p.Results.lognorm;
nbins = p.Results.nbins;
graphics = p.Results.graphics;

x = x(:);   % make sure column vec

if lognorm
    x = log10(x);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[h, e] = histcounts(x, nbins, 'Normalization', 'probability');
e = e(1 : end - 1) + mean(diff(e)) / 2;
% h = movmean(h, 20);
[~, i] = findpeaks(-h);

if lognorm
    sep = 10 ^ mean(e(i));
    % test for normality
    [~, p] = lillietest(x);
    txt = sprintf('thr = %.2f, p(~normality) = %.3f', sep, p);
else
    sep = mean(e(i));
    txt = sprintf('thr = %.2f', sep);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if graphics
    figure
    plot(e, h)
    hold on
    plot([mean(e(i)) mean(e(i))], ylim)
    title(txt)
    ylabel('Norm. Counts')
end


% EOF