function [lRat, iDist, mDist] = cluDist(fet, cluidx, mDist)

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
% 
% TO DO LIST:
%   add option to remove known MU from analysis 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nspk = size(fet, 1);
ncluspk = length(cluidx);
df = size(fet, 2);
if ncluspk < df
    lRat = NaN;
    iDist = NaN;
    mDist = NaN;
    warning('more features then spikes. All params are NaN')
    return
end

% mark spikes which are not cluster members
noiseSpk = setdiff(1:nspk, cluidx);
cluSpk = ismember(1:nspk, cluidx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mahalanobis distances
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 3 
	mDist = mahal(fet, fet(cluidx, :));
end
mCluster = mDist(cluidx); % mahal dist of spikes in cluster
mNoise = mDist(noiseSpk); % mahal dist of all other spikes

% check there is more than one cluster in spike group
if isempty(noiseSpk)
    lRat = 0;
    iDist = Inf;
    warning('only one cluster in group. L-Ratio = 0 and iDist = Inf')
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L ratio
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lRat = sum(1 - chi2cdf(mNoise, df)) / length(cluidx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Isolation distance
% calculate point where mD of other spikes = n of this cell
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ncluspk < nspk / 2  
	[sorted, ~] = sort(mNoise);
	iDist = sorted(ncluspk);
else
    iDist = median(mNoise);
end

end

% EOF