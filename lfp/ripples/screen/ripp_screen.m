function res = ripp_screen(varargin)
% RIPP_SCREEN Compare ripple-detection methods across genotypes.
%
%   res = RIPP_SCREEN(varargin)
%
%   SUMMARY:
%       Runs every method (ripp_methods('screen')) on every session over a short
%       slice, in memory, loading and filtering each session once (methods that
%       share a signal config reuse one prepared signal). Two things come out.
%       First, a compact quality comparison per method: how many events are
%       accepted, their rate, the whitened ripple frequency, and the calibrated
%       threshold. Second, a per-event table of the ACCEPTED NREM events tagged
%       by genotype (mcu_geno: Control, MCU-KO, CAG-MCU-KO), fit with the LME
%       (amp/freq/dur ~ genotype) per method and handed to guiTbl_bar. Nothing is
%       written to a data folder; the results struct is saved once to the
%       manuscript Results directory. Feed res to ripp_screenGui to eyeball where
%       methods disagree.
%
%   INPUTS (Parameter/Value):
%       'basepaths' - <cell> session dirs. {wt_bsl_ripp + mcu_bsl + ra}
%       'win'       - <vec>  slice [start end] (s). {[0 3*3600]}
%       'methods'   - <struct> method configs. {ripp_methods('screen')}
%       'savepath'  - <char> dir for the results .mat. {MCU Results}
%       'flgSave'   - <log>  save res (with the table). {true}
%       'flgLme'    - <log>  fit and print the genotype LMEs. {true}
%       'flgPlot'   - <log>  open guiTbl_bar over the event table. {false}
%       'verbose'   - <log>  per-detection progress. {true}
%
%   OUTPUT:
%       res - <struct> .methods .basepaths .sbjID .genotype .win
%                      .detect{method}(mouse).ripp - all events + .accepted
%                      .pk{mouse,method} - accepted peak times (s)
%                      .meta(mouse,method) - per-run quality
%                      .tbl - per-event table of accepted NREM events.
%
%   DEPENDENCIES:
%       mcu_basepaths, mcu_cfg, basepaths2vars, ripp_methods, ripp_detect,
%       evt_gate, lme_analyse, guiTbl_bar.
%
%   HISTORY:
%       260716 detection-review parameter screen.
%       260719 rebuilt on ripp_detect (shared core); accepted-based; MUA gate.

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addParameter(p, 'basepaths', {}, @iscell);
addParameter(p, 'win', [0 3 * 3600], @isnumeric);
addParameter(p, 'methods', [], @(x) isempty(x) || isstruct(x));
addParameter(p, 'savepath', '', @ischar);
addParameter(p, 'flgSave', true, @islogical);
addParameter(p, 'flgLme', true, @islogical);
addParameter(p, 'flgPlot', false, @islogical);
addParameter(p, 'verbose', true, @islogical);
parse(p, varargin{:});
basepaths = p.Results.basepaths;
win       = p.Results.win;
methods   = p.Results.methods;
savepath  = p.Results.savepath;
flgSave   = p.Results.flgSave;
flgLme    = p.Results.flgLme;
flgPlot   = p.Results.flgPlot;
verbose   = p.Results.verbose;

if isempty(basepaths)
    basepaths = [mcu_basepaths('wt_bsl_ripp'), mcu_basepaths('mcu_bsl'), ...
        mcu_basepaths('ra')];
end
if isempty(methods), methods = ripp_methods('screen'); end
if isempty(savepath)
    savepath = 'D:\OneDrive - Tel-Aviv University\PhD\Slutsky\Manuscripts\MCU\Results';
end
nMice = numel(basepaths);
nM    = numel(methods);
cfg = mcu_cfg;
genoLvl = cfg.lbl.grp;

%% ========================================================================
%  DETECT (mice x methods; each session loaded once, signal reused)
%  ========================================================================

sbjID  = cell(1, nMice);
geno   = cell(1, nMice);
detect = cell(1, nM);
for iMethod = 1:nM
    detect{iMethod} = struct('ripp', cell(1, nMice));
end
pk   = cell(nMice, nM);
meta = repmat(emptyMeta(), nMice, nM);

for iMouse = 1:nMice
    [~, basename] = fileparts(basepaths{iMouse});
    sbjID{iMouse} = strtok(basename, '_');
    geno{iMouse}  = char(mcu_geno(sbjID(iMouse)));
    if verbose
        fprintf('[SCREEN] mouse %d/%d: %s (%s)\n', iMouse, nMice, ...
            sbjID{iMouse}, geno{iMouse});
    end

    v = basepaths2vars('basepaths', basepaths(iMouse), ...
        'vars', {'session', 'spikes', 'spktimes', 'sleep_states', 'units'});
    if ~isfield(v, 'session') || isempty(v.session)
        warning('ripp_screen:noSession', 'no session for %s; skipping', ...
            sbjID{iMouse});
        continue;
    end

    sigCache = struct();
    for iMethod = 1:nM
        met = methods(iMethod);
        key = cfgKey(met);
        sigArg = [];
        if isfield(sigCache, key), sigArg = sigCache.(key); end

        try
            [ripp, aux] = ripp_detect(basepaths{iMouse}, 'met', met, ...
                'win', win, 'v', v, 'sig', sigArg, 'verbose', verbose);
        catch ME
            warning('ripp_screen:detect', '%s / %s failed: %s', ...
                sbjID{iMouse}, met.name, ME.message);
            continue;
        end
        if ~isfield(sigCache, key), sigCache.(key) = aux.sig; end

        % apply the method's QA gate (detect now seeds accepted all-true)
        ripp.accepted = evt_gate(ripp, met.qa);

        % absolute time so ripp_screenGui overlays on the session
        ripp.times    = ripp.times + win(1);
        ripp.peakTime = ripp.peakTime + win(1);

        detect{iMethod}(iMouse).ripp = ripp;
        pk{iMouse, iMethod} = ripp.peakTime(ripp.accepted);
        meta(iMouse, iMethod) = mkMeta(met, ripp, aux, win);
    end
end

res.methods   = methods;
res.basepaths = basepaths;
res.sbjID     = sbjID;
res.genotype  = geno;
res.win       = win;
res.detect    = detect;
res.pk        = pk;
res.meta      = meta;

%% ========================================================================
%  EVENT TABLE (accepted NREM events, genotype-tagged)
%  ========================================================================

amp = []; dur = []; freq = []; freqPeak = []; peakProm = []; gain = [];
methodC = {}; sbjC = {}; genoC = {};
for iMouse = 1:nMice
    for iMethod = 1:nM
        ripp = detect{iMethod}(iMouse).ripp;
        if isempty(ripp) || ~isfield(ripp, 'amp'), continue; end
        keep = ripp.accepted & ripp.state == 'NREM';
        nKeep = sum(keep);
        if nKeep == 0, continue; end
        amp      = [amp;      ripp.amp(keep)];        %#ok<AGROW>
        dur      = [dur;      ripp.dur(keep)];        %#ok<AGROW>
        freq     = [freq;     ripp.freq(keep)];       %#ok<AGROW>
        freqPeak = [freqPeak; ripp.freqPeak(keep)];   %#ok<AGROW>
        peakProm = [peakProm; ripp.peakProm(keep)];   %#ok<AGROW>
        gain     = [gain;     ripp.spkGain(keep)];    %#ok<AGROW>
        methodC  = [methodC;  repmat({methods(iMethod).name}, nKeep, 1)]; %#ok<AGROW>
        sbjC     = [sbjC;     repmat(sbjID(iMouse), nKeep, 1)];  %#ok<AGROW>
        genoC    = [genoC;    repmat(geno(iMouse), nKeep, 1)];   %#ok<AGROW>
    end
end

tbl = table(amp, dur, freq, freqPeak, peakProm, gain, ...
    categorical(methodC), categorical(sbjC), ...
    removecats(categorical(genoC, genoLvl)), ...
    'VariableNames', {'amp', 'dur', 'freq', 'freqPeak', 'peakProm', 'gain', ...
    'method', 'sbjID', 'genotype'});
res.tbl = tbl;

%% ========================================================================
%  REPORT (compact quality per method)
%  ========================================================================

names = {methods.name};
fprintf('\n================ RIPPLE DETECTION SCREEN ================\n');
fprintf('%d mice, window [%g %g] h (accepted NREM events)\n\n', ...
    nMice, win(1) / 3600, win(2) / 3600);

labels = {'accepted (total)', 'rate Hz (mean)', 'freq peak Hz (med)', ...
    'MUA gain (med)', 'thr peak SD (med)', 'chi (med)'};
Q = nan(numel(labels), nM);
for iMethod = 1:nM
    mCol = meta(:, iMethod);
    Q(1, iMethod) = sum([mCol.nAcc], 'omitnan');
    Q(2, iMethod) = mean([mCol.rateHz], 'omitnan');
    Q(3, iMethod) = medianCol(tbl, 'freqPeak', names{iMethod});
    Q(4, iMethod) = medianCol(tbl, 'gain', names{iMethod});
    Q(5, iMethod) = median([mCol.thrPk], 'omitnan');
    Q(6, iMethod) = median([mCol.chi], 'omitnan');
end

fprintf('%-20s', 'metric');
for iMethod = 1:nM, fprintf('%12s', names{iMethod}); end
fprintf('\n');
for iRow = 1:numel(labels)
    fprintf('%-20s', labels{iRow});
    for iMethod = 1:nM, fprintf('%12.3g', Q(iRow, iMethod)); end
    fprintf('\n');
end
fprintf('========================================================\n\n');

%% ========================================================================
%  SAVE (never a data folder)
%  ========================================================================

if flgSave
    if ~isfolder(savepath), savepath = fileparts(mfilename('fullpath')); end
    fname = fullfile(savepath, 'ripp_screen_results.mat');
    save(fname, 'res', '-v7.3');
    if verbose, fprintf('[SCREEN] saved %s\n', fname); end
end

%% ========================================================================
%  GENOTYPE (LME per method + optional bar GUI)
%  ========================================================================

if flgLme && ~isempty(tbl)
    params = {'amp', 'log-normal'; 'freq', 'normal'; 'dur', 'log-normal'};
    for iMethod = 1:nM
        tblM = tbl(tbl.method == names{iMethod}, :);
        fprintf('\n---- genotype LME: method %s ----\n', names{iMethod});
        for iPar = 1:size(params, 1)
            frml = sprintf('%s ~ genotype + (1|sbjID)', params{iPar, 1});
            lme_analyse(tblM, frml, 'dist', params{iPar, 2}, 'verbose', true);
        end
    end
end

if flgPlot && ~isempty(tbl)
    guiTbl_bar(tbl, 'xVar', 'genotype', 'yVar', 'freqPeak', 'grpVar', 'method');
end

end     % EOF


% =========================================================================
%  LOCALS
% =========================================================================
function key = cfgKey(met)
% signal-config identity: methods with the same key share one prepared signal
key = matlab.lang.makeValidName(sprintf('%s_%s_%d_%s', met.chMode, ...
    num2str(met.passband), met.detectMet, met.zMet));
end

% -------------------------------------------------------------------------
function m = mkMeta(met, ripp, aux, win)
% per-run quality summary from one detection
acc = ripp.accepted;
m.name    = met.name;
m.rippCh  = aux.rippCh;
m.nAcc    = sum(acc);
m.rateHz  = sum(acc) / (win(2) - win(1));
m.thrPk   = ripp.info.thr(2);
m.chi     = ripp.info.chi;
end

% -------------------------------------------------------------------------
function m = emptyMeta()
% a meta row for a mouse that was skipped or failed
m = struct('name', '', 'rippCh', [], 'nAcc', NaN, 'rateHz', NaN, ...
    'thrPk', NaN, 'chi', NaN);
end

% -------------------------------------------------------------------------
function y = medianCol(tbl, var, methodName)
% median of a table column for one method, NaN-safe
idx = tbl.method == methodName;
if ~any(idx), y = NaN; return; end
y = median(tbl.(var)(idx), 'omitnan');
end
