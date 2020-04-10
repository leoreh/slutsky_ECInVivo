function isi = ISI(spikes, varargin)

% wrapper for firing rate functions
%
% INPUT
%   spikes      struct (see getSpikes). required input.
%   basepath    path to recording {pwd}
%   ref         refractory period [s] {0.003}.
%   thr         thresholf for contamination [%] {1}.
%   graphics    logical. plot graphics {true} or not (false)
%   saveFig     logical. save figure {1} or not (0)
%
% OUTPUT
%   isi         vector of percent contamination
%               normMethod, normWin
%
% 15 nov 19 LH.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addOptional(p, 'basepath', pwd);
addOptional(p, 'ref', 0.003, @isscalar);
addOptional(p, 'thr', 1, @isscalar);
addOptional(p, 'graphics', 1, @islogical);
addOptional(p, 'saveFig', 0, @islogical);

parse(p, varargin{:})
basepath = p.Results.basepath;
ref = p.Results.ref;
thr = p.Results.thr;
graphics = p.Results.graphics;
saveFig = p.Results.saveFig;

nunits = length(spikes.UID);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% isi
isi = zeros(nunits, 1);
for i = 1 : nunits
    sz(i, 1) = length(spikes.times{i});
    isi(i, 1) = sum(diff(spikes.times{i}) < ref) / sz(i, 1) * 100;
end

% su / mu
su = isi < thr;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if graphics
    f = figure;
    % scatter plot
    subplot(1, 2, 1)
    scatter(sort(isi), [1 : nunits] / nunits * 100, rescale(sz, 20, 100), 'k', 'filled')
    hold on
    axis tight
    ax = gca;
    plot([thr thr], ax.YLim, '--k', 'LineWidth', 2)
    xlabel('ISI contamination [%]')
    ylabel('Unit [%]')
    legend({'SZ = No. Spikes', 'Threshold'}, 'Location', 'southeast')

    
    % distribution
    subplot(1, 2, 2)
    binsize = 0.5;
    bins = [0 : binsize : 5];
    h = histogram(isi(:, 1), bins);
    h.EdgeColor = 'none';
    h.FaceColor = 'k';
    hold on
    axis tight
    ax = gca;
    plot([thr thr], ax.YLim, '--k', 'LineWidth', 2)
    ylabel('No. clusters')
    xlabel('ISI contamination [%]')
    txt = sprintf('ISI contamination within %d ms', ref * 1000);
    suptitle(txt)
    
    if saveFig
        filename = 'ISI distribution';
        savePdf(filename, basepath, f)
    end
    
end

% EOF
