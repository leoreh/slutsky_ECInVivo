function res = ripp_screen(varargin)
% RIPP_SCREEN Compare ripple-detection methods on a session cohort.
%
%   res = RIPP_SCREEN(varargin)
%
%   Runs every method (ripp_screenMethods) on every session over a slice, in
%   memory, and prints one compact comparison: how many events each method
%   finds, their rate, how often an event coincides with a multi-unit spike
%   burst (MUA convergence - a label-free quality proxy), and ripple frequency
%   from both the Hilbert phase and the whitened spectral peak. Nothing is
%   written; the canonical <basename>.ripp.mat is never touched. Feed res to
%   ripp_screenGui to eyeball where methods disagree.
%
%   INPUTS (Parameter/Value):
%       'basepaths' - <cell> session dirs.
%                     {mcu_basepaths('wt_bsl_ripp') + ('mcu_bsl')}
%       'win'       - <vec>  slice [start end] (s). {[0 6*3600]}
%       'methods'   - <struct> method configs. {ripp_screenMethods()}
%       'verbose'   - <log>  per-detection progress. {true}
%
%   OUTPUT:
%       res - <struct> .methods .basepaths .sbjID .win
%                      .detect{method}(mouse).{ripp,rippStates,meta}
%                      .pk{mouse,method} - the method's peak times (s), the
%                                          handle ripp_screenGui reads.
%
%   DEPENDENCIES:
%       mcu_basepaths, basepaths2vars, ripp_screenMethods, ripp_screenDetect.
%
%   HISTORY:
%       260716 detection-review parameter screen.

%% ---- arguments ----
p = inputParser;
addParameter(p, 'basepaths', {}, @iscell);
addParameter(p, 'win', [0 6*3600], @isnumeric);
addParameter(p, 'methods', [], @(x) isempty(x) || isstruct(x));
addParameter(p, 'verbose', true, @islogical);
parse(p, varargin{:});
basepaths = p.Results.basepaths;
win       = p.Results.win;
methods   = p.Results.methods;
verbose   = p.Results.verbose;

if isempty(basepaths)
    basepaths = [mcu_basepaths('wt_bsl_ripp'), mcu_basepaths('mcu_bsl')];
end
if isempty(methods), methods = ripp_screenMethods(); end
nMice = numel(basepaths);
nM    = numel(methods);

%% ---- detect: mice x methods (load each session once) ----
sbjID  = cell(1, nMice);
detect = cell(1, nM);
for k = 1:nM, detect{k} = struct('ripp', {}, 'rippStates', {}, 'meta', {}); end
pk = cell(nMice, nM);

for iMouse = 1:nMice
    [~, basename] = fileparts(basepaths{iMouse});
    sbjID{iMouse} = strtok(basename, '_');
    if verbose
        fprintf('[SCREEN] mouse %d/%d: %s\n', iMouse, nMice, sbjID{iMouse});
    end

    v = basepaths2vars('basepaths', basepaths(iMouse), ...
        'vars', {'session', 'sleep_states', 'spikes'});
    if ~isfield(v, 'session') || isempty(v.session)
        warning('ripp_screen:noSession', 'no session for %s; skipping', sbjID{iMouse});
        continue;
    end

    for k = 1:nM
        try
            [ripp, rippStates, meta] = ripp_screenDetect(basepaths{iMouse}, ...
                methods(k), 'win', win, 'v', v, 'verbose', false);
        catch ME
            warning('ripp_screen:detect', '%s / %s failed: %s', ...
                sbjID{iMouse}, methods(k).name, ME.message);
            ripp = struct('peakTime', [], 'times', zeros(0, 2), ...
                'freqPeak', [], 'freq', [], 'amp', [], 'dur', []);
            rippStates = table();
            meta = struct('nEvents', 0, 'rateHz', NaN, 'muaPos', NaN);
        end
        detect{k}(iMouse).ripp       = ripp;
        detect{k}(iMouse).rippStates = rippStates;
        detect{k}(iMouse).meta       = meta;
        pk{iMouse, k} = ripp.peakTime(:);
    end
end

res.methods = methods; res.basepaths = basepaths; res.sbjID = sbjID;
res.win = win; res.detect = detect; res.pk = pk;

%% ---- report: one compact quality table ----
names = {methods.name};

fprintf('\n================ RIPPLE DETECTION SCREEN ================\n');
fprintf('%d mice, window [%g %g] h\n\n', nMice, win(1)/3600, win(2)/3600);

labels = {'events (total)', 'rate Hz (mean)', 'MUA conv % (med)', ...
    'freq peak Hz', 'freq Hilb Hz'};
fmt = {'%14.0f', '%14.3f', '%14.0f', '%14.0f', '%14.0f'};
Q = zeros(numel(labels), nM);
for k = 1:nM
    d = detect{k};
    Q(1, k) = sum(arrayfun(@(x) x.meta.nEvents, d));
    Q(2, k) = mean(arrayfun(@(x) x.meta.rateHz, d), 'omitnan');
    Q(3, k) = 100 * median(arrayfun(@(x) x.meta.muaPos, d), 'omitnan');
    Q(4, k) = median(arrayfun(@(x) median(x.ripp.freqPeak, 'omitnan'), d), 'omitnan');
    Q(5, k) = median(arrayfun(@(x) median(x.ripp.freq, 'omitnan'), d), 'omitnan');
end

fprintf('%-18s', 'metric');
for k = 1:nM, fprintf('%14s', names{k}); end
fprintf('\n');
for r = 1:numel(labels)
    fprintf('%-18s', labels{r});
    for k = 1:nM, fprintf(fmt{r}, Q(r, k)); end
    fprintf('\n');
end

% per-mouse events, so disagreement (and the hardest session) is visible
fprintf('\n-- events per mouse --\n');
fprintf('%-10s', 'mouse');
for k = 1:nM, fprintf('%14s', names{k}); end
fprintf('\n');
for iMouse = 1:nMice
    fprintf('%-10s', sbjID{iMouse});
    for k = 1:nM
        fprintf('%14d', numel(pk{iMouse, k}));
    end
    fprintf('\n');
end
fprintf('========================================================\n\n');

end     % EOF
