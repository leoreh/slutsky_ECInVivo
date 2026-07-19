function varMap = var_load(varMap, basepath, ctx)

% Fill a varMap of recipes with data, in place, for one session.
%
% Iterates a varMap (a struct: name -> recipe) and writes each recipe's loaded
% value back into its entry as .data / .fs, via var_fetch. An entry that already
% holds non-empty .data is left alone, so a partly-filled map tops up cheaply
% (e.g. across viewer preset switches). An entry whose recipe fails to load is
% dropped, with a warning. The recipe and any other subfields are preserved.
%
% The view half lives elsewhere: a guiMap of panels references these entries by
% name. var_load never touches view.
%
% INPUTS
% - varMap          <struct> name -> recipe (see var_recipe). A recipe may also
%                   be a bare 'file.path' string or an inline value.
% - basepath        <char>(opt) session folder. Default pwd.
% - ctx             <struct>(opt) var_ctx cache to share across calls. Default:
%                   a fresh one at basepath.
%
% OUTPUTS
% - varMap          <struct> the same map, each surviving entry now a struct
%                   holding .data / .fs (fs = NaN for a plain read).
%
% DEPENDENCIES
% - var_fetch, var_ctx.
%
% HISTORY
% - 260719          rewritten as the single loader (recipe grammar + var_fetch;
%                   was the map-fill loop of the earlier var_load).

if nargin < 2 || isempty(basepath), basepath = pwd; end
if nargin < 3 || isempty(ctx), ctx = var_ctx(basepath); end

fldNames = fieldnames(varMap);
drop     = false(1, numel(fldNames));
for iFld = 1 : numel(fldNames)
    entry = varMap.(fldNames{iFld});
    if isstruct(entry) && isfield(entry, 'data') && ~isempty(entry.data)
        continue                                % already filled -> skip
    end
    try
        [data, fs, meta] = var_fetch(entry, ctx);
    catch ME
        warning('var_load:fetch', 'dropping "%s": %s', fldNames{iFld}, ME.message);
        drop(iFld) = true;
        continue
    end
    if ~isstruct(entry), entry = struct(); end  % a bare-string / inline entry
    entry.data   = data;
    entry.fs     = fs;
    entry.labels = meta.labels;                 % channel numbers (bin), else []
    varMap.(fldNames{iFld}) = entry;
end

if any(drop), varMap = rmfield(varMap, fldNames(drop)); end

end

% EOF
