function [sep, h, e] = sepBimodel(varargin)

%  finds the separation point of a bimodel distribution according to the
%  trough in a (log10) histogram. tests if the distribution is lognormal
%  via the lillietest.
% 
% INPUT
%   x           vector to be modled  
%   lognorm     logical. log {1} or linear [0] distribution
%   nbins       scalar {100}. number of bins for histcounts
%   smfactor    scalar {20}. smoothing factor
%   graphics    logical. plot figure {1}.
% 
% OUTPUT
%   sep         trough between distribution
% 
% 22 nov 19 LH
% 30 dec 19 LH  smoothing factor
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
addOptional(p, 'nbins', 100, @isnumeric);
addOptional(p, 'smfactor', 20, @isnumeric);
addOptional(p, 'graphics', true, @islogical);

parse(p, varargin{:})
x = p.Results.x;
lognorm = p.Results.lognorm;
nbins = p.Results.nbins;
smfactor = p.Results.smfactor;
graphics = p.Results.graphics;

x = x(:);   % make sure column vec

if lognorm
    x = log10(x);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% OPTION 1: fit bimodal normal distribution
% from: https://www.mathworks.com/help/stats/examples/fitting-custom-univariate-distributions.html
% pStart = 0.5;
% muStart = quantile(x, [0.25 0.75]);
% sigmaStart = sqrt(var(x) - 0.25 * diff(muStart) .^2);
% start = [pStart muStart sigmaStart sigmaStart];
% lb = [0 -Inf -Inf 0 0];
% ub = [1 Inf Inf Inf Inf];
% 
% pdf_normmixture = @(x,p,mu1,mu2,sigma1,sigma2) ...
%     p*normpdf(x,mu1,sigma1) + (1-p)*normpdf(x,mu2,sigma2);
% paramEsts = mle(x, 'pdf', pdf_normmixture, 'start', start, ...
%     'lower', lb, 'upper', ub);

% OPTION 2: via histogram and findpeaks
[h, e] = histcounts(x, nbins, 'Normalization', 'probability');
e = e(1 : end - 1) + mean(diff(e)) / 2;
h = movmean(h, smfactor);
[p, i] = findpeaks(-h, 'MinPeakDistance', 1);
[~, ii] = min(p); 
i = i(ii);

if lognorm
    sep = 10 ^ mean(e(i));
    % test for normality
    [~, p] = lillietest(x, 'Distribution', 'Normal');
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
    ylabel('Probability')
    set(gca, 'TickLength', [0 0])
    box off
end


% EOF