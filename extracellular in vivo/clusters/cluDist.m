function [L, iDist] = cluDist(fet, cluidx, m)

% compute cluster isolation ditance via:
% (1) L ratio (Schmitzer-Torbert and Redish, 2004)
% (2) Isolation distance (Harris et al., 2001)
%
% INPUT:
%   fet         n (spikes) x d (dimensions) array of feature vectors
%   cluidx      index into fet which lists spikes from the cell whose quality is to be evaluated.
%   m           squared mahalanobis distances, default is to calculate them directly
% 
% OUTPUT:
%   L           L ratio for cluster
%   iDist       isolation distance for cluster
% 
% 03 dec 18 LH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nspk = size(fet, 1);
ncluspk = length(cluidx);
df = size(fet, 2);
if ncluspk < df
    L = nan;
    iDist = NaN;
    warning('more features then spikes. L-Ratio and iDist is NaN')
    return
end

% mark spikes which are not cluster members
noiseSpk = setdiff(1:nspk, cluidx);
cluSpk = ismember(1:nspk, cluidx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mahalanobis distances
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 3 
	m = mahal(fet, fet(cluidx, :));
end
mCluster = m(cluidx); % mahal dist of spikes in cluster
mNoise = m(noiseSpk); % mahal dist of all other spikes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L ratio
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L = sum(1 - chi2cdf(mNoise, df)) / length(cluidx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Isolation distance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate point where mD of other spikes = n of this cell
if ncluspk < nspk / 2  
	[sorted, ~] = sort(mNoise);
	iDist = sorted(ncluspk);
else
    iDist = median(mNoise);
end

end

% EOF