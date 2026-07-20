function fr = spk_rate(spktimes, varargin)
% SPK_RATE Bins spike times into a time-resolved firing rate trace.
%
%   fr = SPK_RATE(spktimes, varargin)
%
%   SUMMARY:
%       Thin wrapper over times2rate that owns the on-disk schema of
%       <basename>.fr.mat. Its only product is the rate trace, where time is
%       a genuine within-unit axis, plus the mean rate over a baseline
%       window. Per-state and per-epoch scalars are NOT computed here - they
%       come from spk_byCond, which labels its rows instead of hiding the
%       condition in an unnamed matrix dimension.
%
%       Smoothing is off by default. Passing smet 'MA' or 'GK' applies a
%       7 bin moving average or gaussian kernel, which at the default 60 s
%       binsize approximates the sliding window of Miyawaki et al., Sci.
%       Rep., 2019.
%
%   INPUTS:
%       spktimes - (Cell)  {nUnits x 1} of spike times [s].
%       varargin - Parameter/Value pairs:
%           'binsize'  - (Num)  Bin size, same units as spktimes. {60}
%           'winCalc'  - (Mat)  [1 x 2] calculation window [s]. {[0 Inf]}
%           'winBL'    - (Mat)  [1 x 2] window for mfr [s]. {[0 Inf]}
%           'smet'     - (Char) Smoothing: 'none' | 'MA' | 'GK'. {'none'}
%           'basepath' - (Char) Recording path. {pwd}
%           'flgSave'  - (Log)  Save <basename>.fr.mat. {true}
%
%   OUTPUTS:
%       fr       - (Struct)
%                    .rate [nUnits x nBins] firing rate [Hz]
%                    .t    [1 x nBins]      bin centres [s]
%                    .mfr  [nUnits x 1]     mean rate across winBL [Hz]
%                    .info parameters
%
%   DEPENDENCIES:
%       times2rate, backup_file.
%
%   HISTORY:
%       26 feb 19 LH  as calc_fr
%       260719        reduced to the rate trace; states, gain, ratio, gini
%                     and fano dropped, strd renamed rate and tstamps t.
%
%   See also: TIMES2RATE, SPK_BYCOND

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 'spktimes', @iscell);
addParameter(p, 'binsize', 60, @isscalar);
addParameter(p, 'winCalc', [0, Inf], @isnumeric);
addParameter(p, 'winBL', [0, Inf], @isnumeric);
addParameter(p, 'smet', 'none', @ischar);
addParameter(p, 'basepath', pwd, @ischar);
addParameter(p, 'flgSave', true, @islogical);

parse(p, spktimes, varargin{:});
binsize  = p.Results.binsize;
winCalc  = p.Results.winCalc;
winBL    = p.Results.winBL;
smet     = p.Results.smet;
basepath = p.Results.basepath;
flgSave  = p.Results.flgSave;

% smoothing width [bins], as in Miyawaki
smFactor = 7;


%% ========================================================================
%  RATE
%  ========================================================================

[fr.rate, ~, tCent] = times2rate(spktimes, 'binsize', binsize, ...
    'winCalc', winCalc, 'c2r', true);
fr.t = tCent(:)';

switch smet
    case 'MA'
        fr.rate = movmean(fr.rate, smFactor, 2);

    case 'GK'
        gk = gausswin(smFactor);
        gk = gk / sum(gk);
        for iUnit = 1 : size(fr.rate, 1)
            fr.rate(iUnit, :) = conv(fr.rate(iUnit, :), gk, 'same');
        end
end

% mean rate across the baseline window. bounds are exclusive, matching the
% legacy definition so migrated and freshly computed files agree
blIdx = fr.t > winBL(1) & fr.t < winBL(2);
fr.mfr = mean(fr.rate(:, blIdx), 2, 'omitnan');


%% ========================================================================
%  FINALIZE
%  ========================================================================

fr.info.binsize = binsize;
fr.info.winCalc = winCalc;
fr.info.winBL   = winBL;
fr.info.smet    = smet;
fr.info.runtime = datetime("now");

if flgSave
    [~, basename] = fileparts(basepath);
    frFile = fullfile(basepath, [basename, '.fr.mat']);
    backup_file(frFile);
    save(frFile, 'fr')
end

end     % EOF
