function [val, meta] = guiPath_src(src, ctx)
% GUIPATH_SRC Resolve a data "address" to its raw value (for guiPath_load).
%
%   [val, meta] = GUIPATH_SRC(src, ctx) turns a terse address into the raw MATLAB
%   value it points at. This is the one place that knows HOW to fetch data, so a
%   panel only has to say WHERE its data lives (see guiPath_panel).
%
%   ADDRESS GRAMMAR (src, a char row like 'PREFIX:VAR.f1.f2'):
%       'ws:VAR[.path]'   base-workspace variable VAR, then a dot-path into it
%                         (the varMap style of io/v2tbl.m). e.g. 'ws:ripp.peakTime'
%       'bin:CH'          raw binary channel CH from <basename>.lfp via binary_load
%                         (nCh / fs / bit2uv taken from <basename>.session.mat).
%       'bin:CH>FILE'     the same, but from the binary FILE given after '>' (for
%                         binaries other than .lfp). e.g. 'bin:5>D:\rec\sess.dat'
%       'sleep_sig[:fld]' the assembled sleep signals, loaded once via ed_sigLoad
%                         (so fs / emg_rms / spec stay consistent); a field of the
%                         sSig struct if given, else the whole sSig struct.
%       'FILE[:VAR.path]' any other prefix is a session file <basename>.FILE.mat,
%                         loaded via basepaths2vars (fuzzy match, auto var name),
%                         then VAR (or the sole variable) and an optional dot-path.
%                         e.g. 'sleep_states:ss.bouts.times', 'ed', 'spikes:spikes.times'
%       (computed adapters - spectrogram, ripple prep - use the 'fn:' scheme and
%        are resolved by guiPath_load, not here.)
%
%   A non-char src is treated as an already-materialized value: a struct with a
%   .data field yields that field, anything else is returned as-is (inline data).
%
%   INPUTS:
%       src  - (Char | numeric | cell | struct) the address or inline value.
%       ctx  - (Struct) shared context from guiPath_ctx: .basepath .basename and a
%              .cache (containers.Map) memoizing loaded files across calls.
%
%   OUTPUTS:
%       val  - the resolved value (vector / matrix / cell / struct).
%       meta - (Struct) .fs   best-effort sampling rate [Hz] (NaN if unknown),
%                       .kind resolver branch taken ('ws'|'bin'|'sleepsig'|'file'|'inline').
%
%   DEPENDENCIES:
%       basepaths2vars, binary_load, ed_sigLoad.
%
%   See also guiPath_load, guiPath_panel, guiPath_presets, v2tbl.
%
%   HISTORY:
%       Created: 05 Jul 2026 - declarative redesign (address resolver).
%       Updated: 05 Jul 2026 - 'bin:CH>FILE' for binaries other than .lfp.

meta = struct('fs', NaN, 'kind', 'inline');

% ---- inline (already-materialized) value ------------------------------------
if ~(ischar(src) || (isstring(src) && isscalar(src)))
    if isstruct(src) && isfield(src, 'data')
        val = src.data;
    else
        val = src;
    end
    return;
end
src = char(src);

% ---- split 'PREFIX:REST' (first colon only) ---------------------------------
ci = find(src == ':', 1);
if isempty(ci)
    prefix = src; rest = '';
else
    prefix = src(1:ci - 1); rest = src(ci + 1:end);
end

switch prefix
    case 'ws'
        meta.kind = 'ws';
        [vname, path] = splitFirst(rest);
        base = evalin('base', vname);
        val  = resolvePath(base, path);

    case 'bin'
        [val, meta] = loadBinCh(rest, ctx);       % returns val + meta.fs

    case 'sleep_sig'
        meta.kind = 'sleepsig';
        ss = getSleepSig(ctx);
        meta.fs = ss.fs;
        if isempty(rest), val = ss.sSig; else, val = resolvePath(ss.sSig, rest); end

    otherwise
        meta.kind = 'file';
        [vname, path] = splitFirst(rest);
        S = getFile(prefix, ctx);                 % struct of the file's variables
        fns = fieldnames(S);
        if ~isempty(vname) && isfield(S, vname)
            base = S.(vname);
        else
            base = S.(fns{1});                    % sole / first variable
            if ~isempty(vname), path = joinPath(vname, path); end
        end
        val = resolvePath(base, path);
        if isstruct(base) && isfield(base, 'fs') && isscalar(base.fs), meta.fs = base.fs; end
end
end

% =========================================================================
%  ADDRESS HELPERS (pure)
% =========================================================================

function [head, tail] = splitFirst(path)
% first dotted token vs the remainder: 'ss.bouts.times' -> 'ss', 'bouts.times'
if isempty(path), head = ''; tail = ''; return; end
di = find(path == '.', 1);
if isempty(di), head = path; tail = ''; else, head = path(1:di - 1); tail = path(di + 1:end); end
end

function p = joinPath(head, tail)
if isempty(tail), p = head; else, p = [head '.' tail]; end
end

function val = resolvePath(base, path)
% walk a dot-path into a struct (the io/v2tbl.m pattern). Empty path -> base.
val = base;
if isempty(path), return; end
parts = strsplit(path, '.');
for i = 1:numel(parts)
    val = val.(parts{i});
end
end

% =========================================================================
%  LOADERS (cached in ctx.cache)
% =========================================================================

function S = getFile(token, ctx)
% load a session file once, as a struct of its variables. Prefer the exact
% <basename>.<token>.mat (matches how the presets loaded ed/ripp/session); fall
% back to basepaths2vars' fuzzy match for compound names (e.g. spikes.cellinfo).
key = ['file:' token];
if isKey(ctx.cache, key), S = ctx.cache(key); return; end
exact = fullfile(ctx.basepath, [ctx.basename, '.', token, '.mat']);
if isfile(exact)
    S = load(exact);
else
    S = basepaths2vars('basepaths', {ctx.basepath}, 'vars', {token}, 'flgPrnt', false);
end
ctx.cache(key) = S;
end

function ss = getSleepSig(ctx)
% assembled sleep signals via ed_sigLoad (canonical loader), cached once. Keeps
% fs / emg_rms / spectrogram exactly as the rest of the pipeline sees them.
if isKey(ctx.cache, 'sleepsig'), ss = ctx.cache('sleepsig'); return; end
[~, ~, ~, fs, specAdapter, sSig] = ed_sigLoad(ctx.basepath, 'basename', ctx.basename);
if (isempty(specAdapter)) && isfield(sSig, 'spec')
    specAdapter = struct('s', sSig.spec, 'freq', sSig.spec_freq(:), 'tstamps', sSig.spec_tstamps(:));
end
ss = struct('sSig', sSig, 'fs', fs, 'spec', specAdapter);
ctx.cache('sleepsig') = ss;
end

function [val, meta] = loadBinCh(chSpec, ctx)
% the requested channel(s) of a binary, as a native [nSamples x nCh] matrix.
% The file is <basename>.lfp by default, or the path after '>' in the spec
% (e.g. '5>D:\rec\sess.dat'); nCh / fs come from <basename>.session.mat.
%
% Returns EVERY requested channel, unaveraged and in the file's native class
% (int16 for an .lfp), so a stack panel gets its columns and the whole channel
% costs a quarter of what a double would. The caller decides what to do with
% the columns: guiPath_load averages them for a 'trace' and keeps them for
% 'traces'. Native means raw ADC counts, not microvolts - fine here, since the
% GUI reads shape, and the ripple / ED presets take physical units from their
% own loaders (fn:ripple.*, fn:edLfp), not from bin:.
meta = struct('fs', NaN, 'kind', 'bin', 'ch', []);
gi = find(chSpec == '>', 1);
if isempty(gi)
    chStr = chSpec; binFile = fullfile(ctx.basepath, [ctx.basename, '.lfp']);
else
    chStr = chSpec(1:gi - 1); binFile = strtrim(chSpec(gi + 1:end));
end
ch = str2num(chStr); %#ok<ST2NM>  (accepts '12' or '[1 2 3]')
if isempty(ch), error('guiPath_src:bin', 'bad binary channel spec "%s"', chSpec); end
session = getSession(ctx);
nCh = session.extracellular.nChannels;
fs  = session.extracellular.srLfp;
val = binary_load(binFile, 'duration', Inf, 'fs', fs, 'nCh', nCh, ...
    'start', 0, 'ch', ch, 'downsample', 1, 'outClass', 'native');
meta.fs = fs;
meta.ch = ch(:)';
end

function session = getSession(ctx)
if isKey(ctx.cache, 'session'), session = ctx.cache('session'); return; end
S = getFile('session', ctx);
if isfield(S, 'session'), session = S.session; else, fns = fieldnames(S); session = S.(fns{1}); end
ctx.cache('session') = session;
end
