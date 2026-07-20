function n = ripp_invalidate(basepath, basename, varargin)
% RIPP_INVALIDATE Delete stale derived ripple products for a session.
%
%   n = RIPP_INVALIDATE(basepath, basename, varargin)
%
%   SUMMARY:
%       Deletes the ACCEPTED-ALIGNED ripple products - spike stats + raster,
%       phase, and optionally the per-bout state table - so they cannot be read
%       stale after the .accepted mask changed. The pipeline stages rebuild them
%       (ripp_curate -> rippStates, ripp_analyze -> spikes / phase). Called by
%       ripp_wrapper on a fresh detection and by ripp_curate when curation
%       changes the mask. Returns the number of files deleted. rippStates is kept
%       by default because a curate save has just rebuilt it; a fresh detection
%       passes flgStates = true to clear it too.
%
%       rippMaps is deliberately NOT here. It holds one row per DETECTED event,
%       row-aligned to ripp, so the mask cannot stale it - readers subset it by
%       .accepted, exactly as they do ripp itself. Only a re-detection changes
%       the event list, and that rewrites the file. Deleting it on a mask change
%       would also destroy the map array the curation GUI reads.
%
%   INPUTS:
%       basepath - <char> session directory.
%       basename - <char> file stem.
%       varargin - Parameter/Value:
%           'flgStates' - <log> also delete rippStates. {false}
%
%   OUTPUT:
%       n - <num> number of product files deleted.
%
%   HISTORY:
%       260719c shared invalidation for the detect -> curate -> analyze split.
%       260720  rippMaps dropped from the list: it became an all-events
%               detection product and is no longer accepted-aligned.

p = inputParser;
addRequired(p, 'basepath', @ischar);
addRequired(p, 'basename', @ischar);
addParameter(p, 'flgStates', false, @islogical);
parse(p, basepath, basename, varargin{:});
flgStates = p.Results.flgStates;

tags = {'rippSpks', 'rippSpkMaps', 'rippSpkLfp'};
if flgStates
    tags = ['rippStates', tags];
end

n = 0;
for iTag = 1:numel(tags)
    f = fullfile(basepath, [basename, '.', tags{iTag}, '.mat']);
    if isfile(f)
        delete(f);
        n = n + 1;
    end
end

end     % EOF
