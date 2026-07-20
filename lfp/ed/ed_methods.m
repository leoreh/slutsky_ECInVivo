function met = ed_methods(preset)
% ED_METHODS Detection + QA configuration (met) for the ED pipeline.
%
%   met = ED_METHODS(preset)
%
%   SUMMARY:
%       The one place that defines an ED "method". The ED twin of ripp_methods:
%       edit a field here, not in the stages.
%
%       Detection and QA are separate concerns. Detection turns the signal into
%       candidates and is deliberately permissive and polarity-blind. QA turns
%       candidates plus their features into an .accepted mask (.qa), which
%       evt_gate applies and ed_curate lets you move per mouse.
%
%       Every value here was measured against 47 discharges curated by hand in
%       three raMCU mice; the record, including what was tried and refuted, is
%       in dev/ed_pipeline_rebuild.md. Read it before moving one.
%
%   FIELDS:
%       .name     - <char> short id, stored as provenance.
%       .chMode   - <char> signal to detect on (ed_sigLoad): 'eeg' = sleep_sig
%                          eeg, present on every scored session; 'ripp' = the
%                          channel ripple detection ran on.
%       .passband - <vec>  detection band [lo hi] (Hz). A discharge is defined
%                          by being sharp, so this sits well above the band the
%                          deflection itself lives in (see ed_detect).
%       .thr      - <num>  candidate threshold on |filt|, in units of the
%                          band's own recording-wide robust scale.
%       .limDur   - <vec>  [min max inter] candidate durations (ms).
%       .qa       - <struct> .ranges: per-metric [lo hi], one field per ed
%                          per-event field. NaN or an absent metric passes.
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
%       260721 rebuilt on curated discharges: 60-150 Hz detection, and the gate
%              reduced to the two criteria that carry it.

if nargin < 1 || isempty(preset), preset = 'default'; end

switch preset
    case 'default'
        met = struct('name', 'default', 'chMode', 'eeg', ...
            'passband', [60 150], 'thr', 8, 'limDur', [4 200 40]);

        % Three criteria, one per way of being wrong, each blind to the
        % others' failure mode:
        %   fastZ  sharp        -> not an ordinary slow LFP deflection
        %   posZ   goes up      -> not a negative step artifact
        %   isoZ   stands alone -> not the biggest peak of a noisy stretch
        % Bounds sit near the 5th percentile of the curated discharges (fastZ
        % p05 18, posZ p05 8, isoZ p05 17), which keeps 43 of 47 of them.
        %
        % This is the set worth REVIEWING, not the answer: it leaves 15-80
        % events on a CAG mouse, which is the ~100 Leore curates down to ~20.
        % Controls land in the single digits and are rejected quickly.
        %
        % emg is computed and is a knob in ed_curate, but is NOT bounded here:
        % once these three are applied it removed nothing measurable, and a
        % bound that does no work will only surprise someone later.
        met.qa.ranges = struct('fastZ', [15 Inf], 'posZ', [5 Inf], ...
            'isoZ', [20 Inf]);

    otherwise
        error('ed_methods:preset', 'unknown preset "%s"', preset);
end

end     % EOF
