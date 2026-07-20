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
%       candidates and is deliberately permissive - a missed discharge is gone
%       for good, an extra candidate costs a row. QA turns candidates plus their
%       features into an .accepted mask (.qa), which evt_gate applies and
%       ed_curate lets you move per mouse.
%
%       Every value here was measured, not assumed; the record is in
%       dev/ed_pipeline_rebuild.md. Read it before moving one.
%
%   FIELDS:
%       .name     - <char> short id, stored as provenance.
%       .chMode   - <char> signal to detect on (ed_sigLoad): 'eeg' = sleep_sig
%                          eeg, present on every scored session; 'ripp' = the
%                          channel ripple detection ran on.
%       .passband - <vec>  discharge band [lo hi] (Hz).
%       .hfBand   - <vec>  supra-physiological band [lo hi] (Hz), behind the
%                          .hfRatio contamination metric.
%       .thr      - <num>  candidate threshold, in recording-robust units.
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

if nargin < 1 || isempty(preset), preset = 'default'; end

switch preset
    case 'default'
        met = struct('name', 'default', 'chMode', 'eeg', ...
            'passband', [10 100], 'hfBand', [150 500], ...
            'thr', 6, 'limDur', [4 200 40]);

        % The gate. Each criterion was ablated against the control-to-epileptic
        % rate contrast and earns its place: dropping ampG costs the most
        % (15.3x -> 9.1x), then ampZ (-> 12.0x), then emg (-> 7.8x). hfRatio
        % moves no rate but halves the survivors on an artifact-dominated
        % session. Together they keep 25 of the 29 discharges a curator
        % confirmed and 6 of the 162 they rejected.
        %
        % emg is loose on purpose. The WAKE/NREM boundary sits near 1, so
        % bounding there deletes the wake discharges the pipeline exists to
        % count - one session fell from 1765 events to 39.
        met.qa.ranges = struct('ampG', [8 Inf], 'ampZ', [7 Inf], ...
            'hfRatio', [-Inf 2], 'emg', [-Inf 2]);

    otherwise
        error('ed_methods:preset', 'unknown preset "%s"', preset);
end

end     % EOF
