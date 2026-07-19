function res = ripp_screen(varargin)
% RIPP_SCREEN Compare ripple-detection methods across genotypes.
%
%   res = RIPP_SCREEN(varargin)
%
%   SUMMARY:
%       Runs every method (ripp_screenMethods) on every session over a short
%       slice, in memory, loading and filtering each session once. Two things
%       come out. First, a compact quality comparison: how many events each
%       method finds, their rate, how often an event coincides with a
%       multi-unit spike burst (MUA convergence, a label-free quality proxy),
%       and the whitened ripple frequency. Second, a per-event table tagged by
%       genotype (Control, MCU-KO, CAG:MCU-KO), fit with the manuscript's LME
%       (amp/freq/dur ~ genotype) so the genotype contrast can be read for each
%       method, and handed to guiTbl_bar for inspection. Nothing is written to
%       a data folder; the canonical <basename>.ripp.mat is never touched. The
%       results struct is saved once to the manuscript Results directory.
%
%       The default methods are the shipping pipeline (current) and the same
%       pipeline with its threshold calibrated to each recording's 1/f noise
%       floor (fooof). Comparing the two shows whether a genotype difference in
%       ripple properties survives holding every recording to a common
%       false-positive rate, which matters across the Intan / TDT batch gap.
%
%   INPUTS (Parameter/Value):
%       'basepaths' - <cell> session dirs. {wt_bsl_ripp + mcu_bsl + ra}
%       'win'       - <vec>  slice [start end] (s). {[0 3*3600]}
%       'methods'   - <struct> method configs. {ripp_screenMethods()}
%       'savepath'  - <char> dir for the results .mat. {MCU Results}
%       'flgSave'   - <log>  save res (with the table). {true}
%       'flgLme'    - <log>  fit and print the genotype LMEs. {true}
%       'flgPlot'   - <log>  open guiTbl_bar over the event table. {false}
%       'verbose'   - <log>  per-detection progress. {true}
%
%   OUTPUT:
%       res - <struct> .methods .basepaths .sbjID .genotype .win
%                      .detect{method}(mouse).{ripp,rippStates}
%                      .pk{mouse,method} - peak times (s), read by ripp_screenGui
%                      .meta(mouse,method) - per-run quality
%                      .tbl - per-event table (amp dur freq freqPeak state
%                             method sbjID genotype).
%
%   DEPENDENCIES:
%       mcu_basepaths, mcu_cfg, basepaths2vars, ripp_screenMethods,
%       ripp_screenDetect, lme_analyse, guiTbl_bar.
%
%   HISTORY:
%       260716 detection-review parameter screen.
%       260717 load once; genotype event table; current vs fooof; LME + bar.

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

cfg = mcu_cfg;
if isempty(basepaths)
    basepaths = [mcu_basepaths('wt_bsl_ripp'), mcu_basepaths('mcu_bsl'), ...
        mcu_basepaths('ra')];
end
if isempty(methods)
    methods = ripp_screenMethods();
end
if isempty(savepath)
    savepath = 'D:\OneDrive - Tel-Aviv University\PhD\Slutsky\Manuscripts\MCU\Results';
end
nMice = numel(basepaths);
nM    = numel(methods);
genoLvl = {'Control', 'MCU-KO', 'CAG:MCU-KO'};

%% ========================================================================
%  DETECT (mice x methods; each session loaded once)
%  ========================================================================

sbjID  = cell(1, nMice);
geno   = cell(1, nMice);
detect = cell(1, nM);
for iMethod = 1:nM
    detect{iMethod} = struct('ripp', cell(1, nMice), 'rippStates', cell(1, nMice));
end
pk   = cell(nMice, nM);
meta = repmat(emptyMeta(), nMice, nM);

for iMouse = 1:nMice
    [~, basename] = fileparts(basepaths{iMouse});
    sbjID{iMouse} = strtok(basename, '_');
    geno{iMouse}  = assignGeno(sbjID{iMouse}, cfg);
    if verbose
        fprintf('[SCREEN] mouse %d/%d: %s (%s)\n', iMouse, nMice, ...
            sbjID{iMouse}, geno{iMouse});
    end

    v = basepaths2vars('basepaths', basepaths(iMouse), ...
        'vars', {'session', 'sleep_states', 'spikes'});
    if ~isfield(v, 'session') || isempty(v.session)
        warning('ripp_screen:noSession', 'no session for %s; skipping', ...
            sbjID{iMouse});
        continue;
    end

    try
        [det, metaM] = ripp_screenDetect(basepaths{iMouse}, methods, ...
            'win', win, 'v', v, 'verbose', verbose);
    catch ME
        warning('ripp_screen:detect', '%s failed: %s', sbjID{iMouse}, ...
            ME.message);
        continue;
    end

    for iMethod = 1:nM
        detect{iMethod}(iMouse).ripp       = det(iMethod).ripp;
        detect{iMethod}(iMouse).rippStates = det(iMethod).rippStates;
        pk{iMouse, iMethod}                = det(iMethod).ripp.peakTime(:);
        meta(iMouse, iMethod)              = metaM(iMethod);
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
%  EVENT TABLE (per event, genotype-tagged)
%  ========================================================================

amp = []; dur = []; freq = []; freqPeak = [];
stateC = {}; methodC = {}; sbjC = {}; genoC = {};
for iMouse = 1:nMice
    for iMethod = 1:nM
        ripp = detect{iMethod}(iMouse).ripp;
        if isempty(ripp) || ~isfield(ripp, 'amp') || isempty(ripp.amp)
            continue;
        end
        nEvt = numel(ripp.amp);
        amp      = [amp;      ripp.amp(:)];        %#ok<AGROW>
        dur      = [dur;      ripp.dur(:)];        %#ok<AGROW>
        freq     = [freq;     ripp.freq(:)];       %#ok<AGROW>
        freqPeak = [freqPeak; ripp.freqPeak(:)];   %#ok<AGROW>
        stateC   = [stateC;  cellstr(string(ripp.state(:)))];        %#ok<AGROW>
        methodC  = [methodC; repmat({methods(iMethod).name}, nEvt, 1)]; %#ok<AGROW>
        sbjC     = [sbjC;    repmat(sbjID(iMouse), nEvt, 1)];        %#ok<AGROW>
        genoC    = [genoC;   repmat(geno(iMouse), nEvt, 1)];         %#ok<AGROW>
    end
end

tbl = table(amp, dur, freq, freqPeak, ...
    categorical(stateC), categorical(methodC), categorical(sbjC), ...
    categorical(genoC, genoLvl), ...
    'VariableNames', {'amp', 'dur', 'freq', 'freqPeak', 'state', ...
    'method', 'sbjID', 'genotype'});
res.tbl = tbl;

%% ========================================================================
%  REPORT (compact quality per method)
%  ========================================================================

names = {methods.name};
fprintf('\n================ RIPPLE DETECTION SCREEN ================\n');
fprintf('%d mice, window [%g %g] h\n\n', nMice, win(1) / 3600, win(2) / 3600);

labels = {'events (total)', 'rate Hz (mean)', 'MUA conv % (med)', ...
    'freq peak Hz (med)', 'thr peak SD (med)', 'chi (med)'};
Q = nan(numel(labels), nM);
for iMethod = 1:nM
    mCol = meta(:, iMethod);
    Q(1, iMethod) = sum([mCol.nEvents], 'omitnan');
    Q(2, iMethod) = mean([mCol.rateHz], 'omitnan');
    Q(3, iMethod) = 100 * median([mCol.muaPos], 'omitnan');
    Q(4, iMethod) = medianCol(tbl, 'freqPeak', names{iMethod});
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
    if ~isfolder(savepath)
        savepath = fileparts(mfilename('fullpath'));
    end
    fname = fullfile(savepath, 'ripp_screen_results.mat');
    save(fname, 'res', '-v7.3');
    if verbose, fprintf('[SCREEN] saved %s\n', fname); end
end

%% ========================================================================
%  GENOTYPE (LME per method + optional bar GUI)
%  ========================================================================

if flgLme
    tblN = tbl(tbl.state == 'NREM', :);
    if isempty(tblN), tblN = tbl; end
    params = {'amp', 'log-normal'; 'freq', 'normal'; 'dur', 'log-normal'};
    for iMethod = 1:nM
        tblM = tblN(tblN.method == names{iMethod}, :);
        fprintf('\n---- genotype LME: method %s ----\n', names{iMethod});
        for iPar = 1:size(params, 1)
            frml = sprintf('%s ~ genotype + (1|sbjID)', params{iPar, 1});
            lme_analyse(tblM, frml, 'dist', params{iPar, 2}, 'verbose', true);
        end
    end
end

if flgPlot
    guiTbl_bar(tbl(tbl.state == 'NREM', :), 'xVar', 'genotype', ...
        'yVar', 'freqPeak', 'grpVar', 'method');
end

end     % EOF


% =========================================================================
%  LOCALS
% =========================================================================
function g = assignGeno(sbjID, cfg)
% genotype from subject id: raMCU = acute viral KO, else the manuscript split
if startsWith(sbjID, 'raMCU')
    g = 'CAG:MCU-KO';
elseif ismember(sbjID, cfg.miceMCU)
    g = 'MCU-KO';
else
    g = 'Control';
end
end

% -------------------------------------------------------------------------
function m = emptyMeta()
% a meta row for a mouse that was skipped or failed
m = struct('name', '', 'rippCh', [], 'fs', NaN, 'win', [NaN NaN], ...
    'nEvents', NaN, 'rateHz', NaN, 'muaPos', NaN, 'thrPk', NaN, 'chi', NaN);
end

% -------------------------------------------------------------------------
function y = medianCol(tbl, var, methodName)
% median of a table column for one method (NREM events), NaN-safe
idx = tbl.method == methodName & tbl.state == 'NREM';
if ~any(idx), idx = tbl.method == methodName; end
y = median(tbl.(var)(idx), 'omitnan');
end
