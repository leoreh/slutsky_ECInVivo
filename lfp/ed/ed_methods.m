function met = ed_methods(preset)
% ED_METHODS Detection + curation configuration (met) for the ED pipeline.
%
%   met = ED_METHODS(preset)
%
%   SUMMARY:
%       The one place that defines an ED "method". The ED twin of ripp_methods:
%       edit a field here, not in the stages.
%
%       Three concerns, kept apart:
%         DETECTION turns the signal into candidates. Permissive and
%           polarity-blind on purpose - a candidate costs a row, a missed
%           discharge is gone for good.
%         .qa is a NOISE FILTER, not a definition. It asks only the two
%           questions no discharge can fail: is it sharp, does it stand alone.
%           Its job is to hand ed_clust a few hundred events instead of a few
%           thousand.
%         .clust groups the survivors by waveform shape, and a human names the
%           groups in ed_curate. That is where "is this a discharge" is
%           decided, which is why no shape criterion lives in .qa.
%
%       There used to be a third gate, posZ, requiring the event to rise above
%       its pre-event baseline. It was dropped: it encoded the shape of the
%       discharges in three mice, polarity is layer-dependent (Maslarova et al.
%       2025 report a sharp negative spike in the dendritic layers and a
%       positive slow wave in the pyramidal layer for the SAME event), and only
%       66% of the curated discharges were positive-going anyway. Shape is now
%       the clustering's business. posZ is still measured and reported.
%
%   FIELDS:
%       .name     - <char> short id, stored as provenance.
%       .passband - <vec>  detection band [lo hi] (Hz). A discharge is defined
%                          by being sharp, so this sits well above the band the
%                          deflection itself lives in (see ed_detect).
%       .thr      - <num>  candidate threshold on |filt|, in units of the
%                          band's own recording-wide robust scale.
%       .limDur   - <vec>  [min max inter] candidate durations (ms).
%       .qa       - <struct> .ranges: per-metric [lo hi], one field per ed
%                          per-event field. NaN or an absent metric passes.
%       .clust    - <struct> ed_clust arguments: .win .nPC .nClust.
%
%   INPUTS:
%       preset - <char> 'default'. {default}
%
%   OUTPUT:
%       met    - <struct> one configuration.
%
%   HISTORY:
%       260720 created with the staged ED pipeline (the met pattern of
%              ripp_methods, replacing ed_wrapper's scattered arguments).
%       260721 rebuilt on curated discharges: 60-150 Hz detection.
%       260721b chMode gone (detection reads one auto-picked .lfp channel);
%              posZ dropped from the gate; .clust added. See
%              dev/ed_pipeline_rebuild.md.

if nargin < 1 || isempty(preset), preset = 'default'; end

switch preset
    case 'default'
        met = struct('name', 'default', 'passband', [60 150], 'thr', 8, ...
            'limDur', [4 200 40]);

        % The noise filter. Two criteria, each blind to the other's failure
        % mode and neither describing a shape:
        %   fastZ  sharp        -> not an ordinary slow LFP deflection
        %   isoZ   stands alone -> not the biggest peak of a noisy stretch
        % Set near the 5th percentile of the curated discharges, so they are
        % permissive by construction: what survives is the set worth SORTING,
        % not the set worth reporting.
        met.qa.ranges = struct('fastZ', [5 Inf], 'isoZ', [5 Inf]);

        % Waveform clustering (ed_clust). Swept against the curated discharges
        % in dev/ed_clustSweep.m and dev/ed_winSweep.m. nClust empty scales the
        % count with the pool (0.65*sqrt(n)) - a fixed count cannot span a pool
        % of 75 and one of 8500, and 12 groups over 8500 leaves every group a
        % mixture. The GUI overrides it per session.
        %
        % .detrend and .norm decide what the components describe. They live
        % here rather than inside ed_clust because a ripple pipeline reusing
        % this curation would want its own answer: L2 divides by the norm over
        % the window, so with events of unequal duration it trades amplitude
        % for length in a way unit peak does not. dev/ed_alignSweep.m could not
        % separate the options on 47 curated discharges - see ed_clust.
        met.clust = struct('win', [-0.05 0.05], 'nPC', 6, 'nClust', [], ...
            'detrend', 'edge', 'norm', 'peak');

    otherwise
        error('ed_methods:preset', 'unknown preset "%s"', preset);
end

end     % EOF
