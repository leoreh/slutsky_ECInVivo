function plotFRdistribution(avgfr, varargin)

% plot the distribution of FR in linear and log scale.
% tests if the distribution is lognormal via the lillietest.
%
% INPUT
%   avgFR       column vector of average FR for each unit
%   saveFig     save figure {1} or not (0)
%   basepath    recording session path {pwd}
%
% 11 jan 19 LH.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addOptional(p, 'saveFig', false, @islogical);
addOptional(p, 'basepath', pwd);

parse(p, varargin{:})
saveFig = p.Results.saveFig;
basepath = p.Results.basepath;

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
xlim([0 10])

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

if saveFig
    filename = 'FR distribution';
    savePdf(filename, basepath, f)
end

end

% EOF