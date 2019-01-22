function f = plotFet(fetMat, varargin)

% plots clusters in feature space
% 
% INPUT
%   fetMat      matrix n x m + 1, where n is the number of spikes and m is
%               the number of features. Last column is the clu ID array
%               includes noise(clu1) and artifact (clu0) spikes. see getFet.
%   clu         clusters to plot
%   fet         vector of features to plot (2 x 1)
%   basepath    recording session path {pwd} to save figure
%   newFig      open new figure or not
%   saveFig     save figure (1) or not {0}
% 
% 12 jan 18 LH.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addOptional(p, 'clu', []);
addOptional(p, 'fet', [1 2], @isnumeric);
addOptional(p, 'basepath', pwd);
addOptional(p, 'newFig', true, @islogical);
addOptional(p, 'saveFig', true, @islogical);

parse(p,varargin{:})
clu = p.Results.clu;
fet = p.Results.fet;
basepath = p.Results.basepath;
newFig = p.Results.newFig;
saveFig = p.Results.saveFig;

if isempty(clu)
    clu = unique(fetMat(:, end));
    clu(clu == 0) = [];     % remove artifact spikes
    clu(clu == 1) = [];     % remove noisy spikes
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if newFig
    f = figure;
else
    f = gcf;
end

hold on
for i = 1 : length(clu)
    idx = find(fetMat(:, end) == clu(i));
    scatter(fetMat(idx, fet(1)), fetMat(idx, fet(2)), '.', 'MarkerEdgeAlpha', 0.05);
    labels{i} = sprintf('clu%d', clu(i));
end
xlim([min(fetMat(:, fet(1))) max(fetMat(:, fet(1)))])
ylim([min(fetMat(:, fet(2))) max(fetMat(:, fet(2)))])
axis tight
xlabel(sprintf('PC%d', fet(1)))
ylabel(sprintf('PC%d', fet(2)))
% legend(labels)
filename = sprintf('fet %s of clu %s', num2str(fet), num2str(clu));
title(filename)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if saveFig
    savePdf(filename, basepath, f)
end

end

% EOF