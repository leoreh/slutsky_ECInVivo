function plotFRdistribution(avgfr, varargin)

% plot the distribution of FR in linear and normal scale.
% tests if the distribution is lognormal via the lillietest.
%
% INPUT
%   avgFR       column vector of FR for each unit (row)
%
% 11 jan 19 LH.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
validate_win = @(win) assert(isnumeric(win) && length(win) == 2,...
    'time window must be in the format [start end]');

p = inputParser;
addOptional(p, 'win', [1 Inf], validate_win);
addOptional(p, 'method', 'avg', @ischar);

parse(p,varargin{:})
win = p.Results.win;
method = p.Results.method;

if win(1) == 0; win(1) = 1; end
if win(2) == Inf; win(2) = size(spkcount.strd, 2); end

nunits = size(spkcount.strd, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = figure;

% linear
subplot(1, 2, 1)
h = histogram(avgfr, 20);
h.EdgeColor = 'none';
h.FaceColor = 'k';
box off
axis tight
xlabel('Firing Rate [Hz]')
ylabel('Number of Units')

% log
subplot(1, 2, 2)
h = histogram(log10(avgfr), 20);
h.EdgeColor = 'none';
h.FaceColor = 'k';
box off
axis tight
xlabel('Firing Rate [log(Hz)]')
ylabel('Number of Units')

% test for normality
[~, p] = lillietest(log10(avgfr));
txt = sprintf('p(~normality) = %.2f', p);
title(txt)

end

% EOF