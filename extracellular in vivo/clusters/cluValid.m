function spikes = cluValid(basepath, spikes, npca, graphics, saveFig)

% wrapper for inspecting and validating clusters
% computes:
%   (1) cluster isolation ditance via:
%       (a) L ratio (Schmitzer-Torbert and Redish, 2004)
%       (b) Isolation distance (Harris et al., 2001)
%   (2) ISI contamination
%
% INPUT:
%   basepath        path to recording
%   spikes          struct (see getSpikes)
%   npca            no. PCA features for mahal dist {12}.
%                   If sorted with neurosuite should be #spk group * 3
%   graphics        logical. plot graphics {true} or not (false)
%   saveFig         save figure {1} or not (0)
%
% OUTPUT:           additional fields to struct (all 1 x nunits):
%   su              SU (1) or MU (0)
%   L                         
%   iDist            
%   ISIratio                 
%
% CALLS:
%   getSpikes       optional (if not given as input)
%   getFet          PCA feature array
%   cluDist         L ratio and iDist
%   plotCluster     plot waveform and ISI histogram
%
% 03 dec 18 LH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nargs = nargin;
if nargs < 1 || isempty(basepath)
    basepath = pwd;
end
if nargs < 2 || isempty(spikes)
    warning('spikes will be loaded from %s', basepath)
    spikes = getSpikes(basepath);
end
if nargs < 3 || isempty(npca)
    npca = 12;
end
if nargs < 4 || isempty(graphics)
    graphics = 1;
end
if nargs < 5 || isempty(saveFig)
    saveFig = true;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate cluster separation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fet = getFet(basepath);

for i = 1 : length(fet)
    nclu = unique(fet{i}(:, end));
    nclu = length(nclu(nclu > 1));
    L{i} = zeros(nclu, 1);
    iDist{i} = zeros(nclu ,1);
    for j = 1 : nclu
        cluidx = find(fet{i}(:, end) == spikes.cluID(j));
        [L{i}(j), iDist{i}(j)] = cluDist(fet{i}(:, 1 : npca), cluidx);
    end
end

% concatenate cell
spikes.L = cat(1, L{:});
spikes.iDist = cat(1, iDist{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate cluster separation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
spikes.ISIratio = zeros(length(spikes.UID), 1);
for i = 1 : length(spikes.UID)
    spikes.ISIratio(i) = sum(diff(spikes.times{i}) < 0.002) / length(spikes.times{i}) * 100;
end

% for i = 1 : length(spikes.UID)
%     spikes.ISIratio(i) = histcounts(diff(spikes.times{i}), [0 0.002]) /...
%         histcounts(diff(spikes.times{i}), [0 0.02]) * 100;
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% estimate SU or MU
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ISIrmv = sum((spikes.ISIratio < 0.1));
Lrmv = sum((spikes.L < 0.05));
iDrmv = sum((spikes.iDist > 20));
iDLrmv = sum((spikes.iDist > 20 & spikes.L < 0.05));
LniD = sum((spikes.iDist > 20 & spikes.L < 0.05));
ISILiD = sum(spikes.L < 0.05 & spikes.iDist > 20 & spikes.ISIratio < 0.1);
ISIid = sum(spikes.iDist > 20 & spikes.ISIratio < 0.1);

spikes.su = spikes.L < 0.05 & spikes.iDist > 20 & spikes.ISIratio < 0.1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if graphics
    
    plotCluster(basepath, spikes, [], saveFig)
    
   
end

end

% EOF