function spikes = cluVal(basepath, spikes, npca, graphics, saveFig, saveVar)

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
%   saveVar         save variables (update spikes and save su)
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
    spikes = getSpikes('basepath', basepath);
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
if nargs < 6 || isempty(saveVar)
    saveVar = true;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate cluster separation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fet = getFet(basepath);

if all(isfield(spikes, {'mDist', 'iDist', 'L'}))
    warning('separation matrices already calculated. Skipping cluDist');
else
    for i = 1 : length(fet)
        
        % disgard noise and artifact spikes
        idx = find(fet{i}(:, end) <= 1);
        fet{i}(idx, :) = [];
        
        % calculate separation matrices
        clu = unique(fet{i}(:, end));
        clu = clu(clu > 1);
        L{i} = zeros(length(clu), 1);
        iDist{i} = zeros(length(clu) ,1);
        for j = 1 : length(clu)
            cluidx = find(fet{i}(:, end) == clu(j));
            [L{i}(j, 1), iDist{i}(j, 1), mDist{i}{j}] = cluDist(fet{i}(:, 1 : npca), cluidx);
        end
    end
    
    % concatenate cell
    spikes.L = cat(1, L{:});
    spikes.iDist = cat(1, iDist{:});
    spikes.mDist = mDist;  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate cluster separation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
spikes.ISIratio = zeros(length(spikes.UID), 1);
for i = 1 : length(spikes.UID)
    spikes.ISIratio(i, 1) = sum(diff(spikes.times{i}) < 0.003) / length(spikes.times{i}) * 100;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% estimate SU or MU
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
spikes.su = spikes.iDist > 20 & spikes.ISIratio < 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save spikes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if saveVar
    [~, filename] = fileparts(basepath);
    save([filename '.spikes.cellinfo.mat'], 'spikes')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if graphics  
    plotCluster(basepath, spikes, [], saveFig)     
end

end

% EOF