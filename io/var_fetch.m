function [data, fs, meta] = var_fetch(recipe, ctx)

% Resolve one recipe to its raw value: fetch -> transform chain -> path walk.
%
% The single-recipe primitive under var_load. It is view-blind: it returns the
% raw value and its sampling rate, never a shaped panel. Reuse consumers (e.g.
% as_prepSig, ripp_wrapper) can call it directly, one signal at a time.
%
% A recipe is built by var_recipe. Two shorthands are accepted: a bare string
% 'file.path' is read as a matvar; any other non-recipe value is taken inline.
%
% INPUTS
% - recipe          <struct | char> a var_recipe, a 'file.path' string, or an
%                   inline value.
% - ctx             <struct>(opt) var_ctx cache. Default: a fresh one at pwd.
%
% OUTPUTS
% - data            the raw value (vector / matrix / cell / struct).
% - fs              <num> sampling rate [Hz]; NaN when the source carries none.
% - meta            <struct> .kind and .labels (channel numbers, for a bin).
%
% DEPENDENCIES
% - binary_load, basepaths2vars (fetch); processEMG, calc_spec, ripp_sigPrep,
%   iosr.dsp.sincFilter (transforms).
%
% HISTORY
% - 260719          created (unified var_* I/O layer; absorbs guiPath_src and
%                   the loading half of guiPath_load).

if nargin < 2 || isempty(ctx), ctx = var_ctx(pwd); end
recipe = toRecipe(recipe);


%% ========================================================================
%  FETCH
%  ========================================================================

meta = struct('kind', recipe.kind, 'labels', []);
switch recipe.kind
    case 'matvar'
        [data, fs] = fetchMatvar(recipe, ctx);
    case 'matfield'
        [data, fs] = fetchMatfield(recipe, ctx);
    case 'bin'
        [data, fs, meta.labels] = fetchBin(recipe, ctx);
    case 'ws'
        [data, fs] = fetchWs(recipe);
    case 'value'
        data = recipe.data;
        fs   = NaN;
    otherwise
        error('var_fetch:kind', 'unknown recipe kind "%s"', recipe.kind);
end


%% ========================================================================
%  TRANSFORM CHAIN, PATH, FS OVERRIDE
%  ========================================================================
% Each transform is its own adapter (a switch, not a dynamic feval, since the
% ops have different signatures). Add a transform -> add a case in applyTransform.

for iStep = 1 : numel(recipe.transform)
    st = recipe.transform(iStep);
    [data, fs] = applyTransform(st.op, data, fs, st.args);
end

if ~isempty(recipe.path)
    data = walkPath(data, recipe.path);
end

if ~isempty(recipe.fs)                          % a read that names its own rate
    fs = recipe.fs;
end

end


% =========================================================================
%  FETCHERS
% =========================================================================

function [data, fs] = fetchMatvar(recipe, ctx)
% the wrapper variable of <basename>.file.mat (named var, else sole/first). The
% dot-path is walked once, later, by var_fetch (after any transform).
S = loadMat(recipe.file, ctx);
if ~isempty(recipe.var) && isfield(S, recipe.var)
    data = S.(recipe.var);
else
    fns  = fieldnames(S);
    data = S.(fns{1});                          % sole / first variable
end
fs = NaN;                                       % signals come from bin/matfield
end


function [data, fs] = fetchMatfield(recipe, ctx)
% a named top-level field of a -struct .mat (e.g. sleep_sig eeg/emg/emg_rms).
% A cellstr field packs those fields into one struct (the spec adapter's shape).
S   = loadMat(recipe.file, ctx);
fld = recipe.field;
if iscell(fld)
    data = struct();
    for iFld = 1 : numel(fld)
        data.(fld{iFld}) = S.(fld{iFld});
    end
else
    data = S.(fld);
end
fs = NaN;
if isfield(S, 'fs') && isscalar(S.fs), fs = S.fs; end
end


function [data, fs, labels] = fetchBin(recipe, ctx)
% channel(s) of a binary via binary_load. nCh / fs from <basename>.session.mat;
% 'native' keeps raw int16 (for a stack), 'double' scales by bit2uv. Several
% channels average into one trace when 'average' is set. ch passes through
% unchanged (binary_load is 1-based).
session = getSession(ctx);
fs  = session.extracellular.srLfp;
nCh = session.extracellular.nChannels;
if isAbsPath(recipe.file)
    binFile = recipe.file;
else
    binFile = fullfile(ctx.basepath, [ctx.basename, '.', recipe.file]);
end
ch     = recipe.ch;
labels = ch;
key = sprintf('bin:%s|%s|%s', binFile, mat2str(ch), recipe.outClass);
if isKey(ctx.cache, key)
    data = ctx.cache(key);
else
    args = {'fs', fs, 'nCh', nCh, 'duration', Inf, 'start', 0, ...
        'ch', ch, 'downsample', 1, 'outClass', recipe.outClass};
    if ~isempty(recipe.bit2uv), args = [args, {'bit2uv', recipe.bit2uv}]; end
    data = binary_load(binFile, args{:});
    ctx.cache(key) = data;
end
if recipe.average && size(data, 2) > 1
    data = mean(data, 2);
end
end


function [data, fs] = fetchWs(recipe)
% a base-workspace variable (the GUI Load-from-workspace source); the dot-path
% is walked later by var_fetch
data = evalin('base', recipe.var);
fs   = NaN;
end


% =========================================================================
%  TRANSFORMS (each op is its own adapter)
% =========================================================================

function [data, fs] = applyTransform(op, data, fs, args)
switch op
    case 'eegSub'
        [data, fs] = tfEegSub(data, fs, args);      % as_prepSig eeg path
    case 'emg'
        [data, fs] = tfEmg(data, fs, args);         % as_prepSig emg path
    case 'emgRms'
        data = processEMG(data(:), fs, 1);
        fs   = 1;
    case 'spec'
        data = calc_spec('sig', double(data(:)), 'fs', fs, 'graphics', false, ...
            'saveVar', false, 'force', true, args{:});
        fs = NaN;
    case 'rippPrep'
        [data, fs] = tfRippPrep(data, fs, args);
    otherwise
        error('var_fetch:transform', 'unknown transform "%s"', op);
end
end


function [sig, fsTarget] = tfEegSub(sig, fsIn, args)
% low-pass (iosr sinc) then subsample, matching as_prepSig's eeg branch
[fsTarget, cf] = deal(args{:});
sig = double(sig(:));
if fsIn ~= fsTarget
    if ~isempty(cf), sig = iosr.dsp.sincFilter(sig, cf / (fsIn / 2)); end
    r   = fsIn / fsTarget;
    sig = sig(r : r : numel(sig));
end
sig = sig(:);
end


function [emg, fsTarget] = tfEmg(sig, fsIn, args)
% spline-resample to fsTarget then band-filter, matching as_prepSig's emg branch
[fsTarget, cf] = deal(args{:});
tRaw = (1 : numel(sig)) / fsIn;
tNew = (1 : floor(tRaw(end) * fsTarget)) / fsTarget;
emg  = interp1(tRaw, double(sig(:)), tNew, 'spline')';
if ~isempty(cf), emg = iosr.dsp.sincFilter(emg, cf / (fsTarget / 2)); end
emg  = emg(:);
end


function [rs, fs] = tfRippPrep(sig, fs, args)
% ripp_sigPrep on a filtered/averaged channel (the ripple curation signal)
passband = args{1};
zMet = 'adaptive';
if numel(args) >= 2 && ~isempty(args{2}), zMet = args{2}; end
rs = ripp_sigPrep(double(sig(:)), fs, 'passband', passband, 'zMet', zMet);
end


% =========================================================================
%  SMALL HELPERS (pure)
% =========================================================================

function recipe = toRecipe(recipe)
% accept a var_recipe struct, a bare 'file.path' string, or an inline value
if ischar(recipe) || (isstring(recipe) && isscalar(recipe))
    recipe = strRecipe(char(recipe));
elseif ~isstruct(recipe) || ~isfield(recipe, 'kind')
    recipe = var_recipe('value', 'data', recipe);
end
end


function recipe = strRecipe(str)
% split a bare-string read 'file.path' into a matvar recipe
di = find(str == '.', 1);
if isempty(di)
    recipe = var_recipe('matvar', 'file', str);
else
    recipe = var_recipe('matvar', 'file', str(1 : di - 1), ...
        'path', str(di + 1 : end));
end
end


function S = loadMat(fileToken, ctx)
% load <basename>.token.mat once (exact name, else a fuzzy basepaths2vars match)
key = ['mat:', fileToken];
if isKey(ctx.cache, key), S = ctx.cache(key); return; end
exact = fullfile(ctx.basepath, [ctx.basename, '.', fileToken, '.mat']);
if isfile(exact)
    S = load(exact);
else
    S = basepaths2vars('basepaths', {ctx.basepath}, 'vars', {fileToken}, ...
        'flgPrnt', false);
end
ctx.cache(key) = S;
end


function session = getSession(ctx)
S = loadMat('session', ctx);
if isfield(S, 'session')
    session = S.session;
else
    fns = fieldnames(S);
    session = S.(fns{1});
end
end


function tf = isAbsPath(f)
% a bin file that is an absolute path (a foreign binary) vs a session token
tf = ~isempty(regexp(f, '^[A-Za-z]:', 'once')) || strncmp(f, '\\', 2) ...
    || strncmp(f, '/', 1);
end


function val = walkPath(base, path)
% walk a dot-path into a struct; empty path -> base
val = base;
if isempty(path), return; end
parts = strsplit(path, '.');
for iPart = 1 : numel(parts)
    val = val.(parts{iPart});
end
end

% EOF
