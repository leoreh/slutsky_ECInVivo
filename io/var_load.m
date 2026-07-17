function varMap = var_load(varMap, basepath)

% Fill a varMap of recipes with data for one session.
%
% Reads each entry's recipe and writes the loaded value back into that entry
% (.data, .fs), leaving the recipe and any .view subfield in place. An entry
% that already holds data is skipped, so a partly-filled map tops up cheaply.
% Every .mat is read once per call.
%
% INPUTS
% - varMap          <struct> name -> recipe (or an already-filled entry).
%                   A recipe is a read or a binary, told apart by 'args':
%                       'fr.mfr'   a field of <basename>.fr.mat.
%                       struct('file','lfp', 'args',{{'ch',1}})   a binary
%                       channel; args pass to binary_load (session gives
%                       default nCh / fs, args override). Add 'fn' to transform
%                       the channel ('spec' | 'emg' | 'emgRms'), 'fnArgs' its
%                       parameters.
% - basepath        <char>(opt) session folder. Default pwd.
%
% OUTPUTS
% - varMap          <struct> the same map, each entry now holding .data / .fs.
%                   fs is NaN for a read, the sampling rate for a signal.
%
% DEPENDENCIES
% - binary_load, calc_spec, processEMG, iosr.dsp.sincFilter (emg transforms).
%
% HISTORY
% - 260717          created (unified var_* I/O layer).

if nargin < 2 || isempty(basepath), basepath = pwd; end
[~, basename] = fileparts(basepath);

cache    = containers.Map;                  % read each .mat once per call
fldNames = fieldnames(varMap);
for iFld = 1 : numel(fldNames)
    entry = varMap.(fldNames{iFld});
    if isstruct(entry) && isfield(entry, 'data')
        continue                            % already filled -> skip
    end
    if ischar(entry) || isstring(entry)
        entry = strToRecipe(char(entry));   % bare-string read -> flat recipe
    end

    if isfield(entry, 'args')               % binary (binary_load)
        [sig, fs] = loadBin(entry, basepath, basename, cache);
        if isfield(entry, 'fn')
            fnArgs = {};
            if isfield(entry, 'fnArgs'), fnArgs = entry.fnArgs; end
            [entry.data, entry.fs] = applyFn(entry.fn, sig, fs, fnArgs);
        else
            entry.data = sig;
            entry.fs   = fs;
        end
    else                                    % read a .mat field
        entry.data = readField(entry, basepath, basename, cache);
        entry.fs   = NaN;
    end

    varMap.(fldNames{iFld}) = entry;
end

end


% =========================================================================
%  HELPERS
% =========================================================================

function recipe = strToRecipe(str)
% split a bare-string read 'file.path' into a flat recipe
di = find(str == '.', 1);
if isempty(di)
    recipe = struct('file', str, 'path', '');
else
    recipe = struct('file', str(1 : di - 1), 'path', str(di + 1 : end));
end
end


function val = readField(recipe, basepath, basename, cache)
% a field of <basename>.file.mat: the variable inside it, then the dot-path
matStruct = loadMat(recipe.file, basepath, basename, cache);
val = mainVar(matStruct, recipe.file);
if isfield(recipe, 'path') && ~isempty(recipe.path)
    parts = strsplit(recipe.path, '.');
    for iPart = 1 : numel(parts)
        val = val.(parts{iPart});
    end
end
end


function [sig, fs] = loadBin(recipe, basepath, basename, cache)
% channel(s) of <basename>.file via binary_load. The session gives the default
% nCh / fs; recipe.args override, so a foreign binary (emg.dat) loads by passing
% its own geometry. Several channels average into one trace.
session = getSession(basepath, basename, cache);
opt = struct('fs', session.extracellular.srLfp, ...
    'nCh', session.extracellular.nChannels, 'duration', Inf, 'downsample', 1);
opt = applyArgs(opt, recipe.args);
fname = fullfile(basepath, [basename, '.', recipe.file]);
nv  = nvPairs(opt);
sig = binary_load(fname, nv{:});
if size(sig, 2) > 1
    sig = mean(sig, 2);
end
fs = opt.fs;
end


function [data, fs] = applyFn(fn, sig, fsIn, fnArgs)
% run a transform on a loaded channel, passing fnArgs. A switch, not a dynamic
% feval: each transform has its own signature, so the switch is the adapter.
% Add a transform -> add a case.
switch fn
    case 'spec'
        data = calc_spec('sig', sig, 'fs', fsIn, 'graphics', false, ...
            'saveVar', false, 'force', true, fnArgs{:});
        fs = NaN;

    case 'emg'
        [data, fs] = emgResample(sig, fsIn, fnArgs);

    case 'emgRms'
        [emg, fsEmg] = emgResample(sig, fsIn, fnArgs);
        data = processEMG(emg, fsEmg, 1);
        fs = 1;

    otherwise
        error('var_load:fn', 'unknown transform "%s"', fn);
end
end


function [emg, fsTarget] = emgResample(sig, fsIn, fnArgs)
% resample a raw EMG channel to fsTarget and band-filter it, matching
% as_prepSig (spline to the target length, then the iosr sinc filter)
p = applyArgs(struct('fsTarget', 1250, 'cf', [80 450]), fnArgs);
fsTarget = p.fsTarget;
tRaw = (1 : numel(sig)) / fsIn;
tNew = (1 : floor(tRaw(end) * fsTarget)) / fsTarget;
emg  = interp1(tRaw, double(sig(:)), tNew, 'spline')';
if ~isempty(p.cf)
    emg = iosr.dsp.sincFilter(emg, p.cf / (fsTarget / 2));
end
emg = emg(:);
end


function matStruct = loadMat(token, basepath, basename, cache)
% load <basename>.token.mat once, exact name first then a wildcard
if isKey(cache, token)
    matStruct = cache(token);
    return
end
fname = fullfile(basepath, [basename, '.', token, '.mat']);
if ~isfile(fname)
    hits = dir(fullfile(basepath, ['*', token, '*.mat']));
    if isempty(hits)
        error('var_load:file', 'no "%s" file in %s', token, basepath);
    end
    fname = fullfile(basepath, hits(1).name);
end
matStruct = load(fname);
cache(token) = matStruct;
end


function val = mainVar(matStruct, token)
% the variable inside a .mat: the one named token, else the first
if isfield(matStruct, token)
    val = matStruct.(token);
else
    fldNames = fieldnames(matStruct);
    val = matStruct.(fldNames{1});
end
end


function session = getSession(basepath, basename, cache)
session = mainVar(loadMat('session', basepath, basename, cache), 'session');
end


function opt = applyArgs(opt, args)
% fold a name-value cell into the option struct (args win over defaults)
for iArg = 1 : 2 : numel(args) - 1
    opt.(args{iArg}) = args{iArg + 1};
end
end


function nv = nvPairs(s)
% a struct to a name-value cell for a varargin call
fldNames = fieldnames(s);
nv = cell(1, 2 * numel(fldNames));
nv(1 : 2 : end) = fldNames;
nv(2 : 2 : end) = struct2cell(s);
end

% EOF
