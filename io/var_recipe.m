function recipe = var_recipe(kind, varargin)

% Build one validated data recipe for var_load / var_fetch.
%
% A recipe says WHERE a signal lives and, optionally, how to transform it -
% never how to draw it. var_load fills its .data / .fs in place; var_fetch
% resolves one on its own. A view (a guiMap panel) references a recipe by name;
% it does not carry one.
%
% Resolution per recipe is: fetch (by kind) -> transform chain -> path walk.
% The loader stays view-blind: it returns the raw value, not a shaped panel.
%
% EXAMPLES
% - var_recipe('matvar', 'file','sleep_states', 'var','ss', 'path','bouts.times')
%   a field of <basename>.sleep_states.mat (wrapper variable ss, then a dot-path).
% - var_recipe('matfield', 'file','sleep_sig', 'field','emg_rms', 'fs',1)
%   a named top-level field of a -struct .mat (sleep_sig has no wrapper var).
% - var_recipe('bin', 'file','lfp', 'ch',[5 6 7 8], 'average',true)
%   channels 5-8 of <basename>.lfp, averaged to one trace.
% - var_recipe('bin', 'file','lfp', 'ch',rippCh, 'bit2uv',0.195, ...
%       'transform',{'bandpass',{[80 250]}})
%   the ripple channel, band-passed to the ripple band ('rippPrep' instead runs
%   the full detection prep and returns its struct, e.g. with 'path','z').
%
% INPUTS
% - kind            <char> 'matvar' | 'matfield' | 'bin' | 'ws' | 'value'.
%
% NAME-VALUE (used per kind; irrelevant ones are ignored)
% - file            <char> matvar/matfield: file token <basename>.file.mat.
%                          bin: 'lfp' (default) or an absolute binary path.
% - var             <char> matvar: wrapper variable inside the .mat.
% - field           <char | cellstr> matfield: top-level field(s); a cellstr
%                          packs those fields into a struct.
% - ch              <num>  bin: channel(s), passed to binary_load unchanged
%                          (1-based, as binary_load expects).
% - average         <logical> bin: mean across channels. Default true.
% - outClass        <char> bin: 'double' (scaled by bit2uv) | 'native'
%                          (raw int16). Default 'double'.
% - bit2uv          <num>  bin: scaling for binary_load. Default [] (its rule).
% - data            <any>  value: an already-materialized value.
% - path            <char> dot-path walked on the final value. Default ''.
% - transform       <cell> op/args pairs run after the fetch, e.g. {'spec',{}}
%                          or {'eegSub',{eegFs,fs,cf}, 'spec',{}}. Default {}.
% - fs              <num>  fs to report for a read that carries none. Default [].
%
% OUTPUTS
% - recipe          <struct> .kind + the kind's fields + .path .transform .fs.
%
% DEPENDENCIES
% - none (a pure constructor). See var_fetch for what consumes a recipe.
%
% HISTORY
% - 260719          created (unified var_* I/O layer).


%% ========================================================================
%  ARGUMENTS
%  ========================================================================

kind = lower(char(kind));

p = inputParser;
p.FunctionName = 'var_recipe';
addParameter(p, 'file',      '');
addParameter(p, 'var',       '');
addParameter(p, 'field',     {});
addParameter(p, 'ch',        []);
addParameter(p, 'average',   true,     @islogical);
addParameter(p, 'outClass',  'double', @ischar);
addParameter(p, 'bit2uv',    []);
addParameter(p, 'data',      []);
addParameter(p, 'path',      '',       @ischar);
addParameter(p, 'transform', {},       @iscell);
addParameter(p, 'fs',        []);
parse(p, varargin{:});
a = p.Results;


%% ========================================================================
%  BUILD
%  ========================================================================
% Build field-by-field (never struct(name, value, ...): a cell / array value
% would spawn a struct array).

recipe = struct('kind', kind);
switch kind
    case 'matvar'
        assert(~isempty(a.file), 'var_recipe:matvar', 'matvar needs a file');
        recipe.file = a.file;
        recipe.var  = a.var;

    case 'matfield'
        assert(~isempty(a.file) && ~isempty(a.field), ...
            'var_recipe:matfield', 'matfield needs a file and a field');
        recipe.file  = a.file;
        recipe.field = a.field;

    case 'bin'
        assert(~isempty(a.file) && ~isempty(a.ch), ...
            'var_recipe:bin', 'bin needs a file and ch');
        assert(all(a.ch >= 1), 'var_recipe:ch', ...
            'bin ch is 1-based (binary_load requires ch >= 1)');
        recipe.file     = a.file;
        recipe.ch       = a.ch(:)';
        recipe.average  = a.average;
        recipe.outClass = lower(a.outClass);
        recipe.bit2uv   = a.bit2uv;

    case 'ws'
        assert(~isempty(a.var), 'var_recipe:ws', 'ws needs a var');
        recipe.var = a.var;

    case 'value'
        recipe.data = a.data;

    otherwise
        error('var_recipe:kind', 'unknown kind "%s"', kind);
end

recipe.path      = a.path;
recipe.transform = normTransform(a.transform);
recipe.fs        = a.fs;

end


% =========================================================================
%  HELPERS
% =========================================================================

function t = normTransform(c)
% a {'op1',{args1}, 'op2',{args2}, ...} cell -> a 1xN struct array .op .args
if isempty(c)
    t = struct('op', {}, 'args', {});
    return
end
assert(mod(numel(c), 2) == 0, 'var_recipe:transform', ...
    'transform must be op/args pairs');
n = numel(c) / 2;
t = struct('op', cell(1, n), 'args', cell(1, n));
for iOp = 1 : n
    t(iOp).op   = c{2 * iOp - 1};
    t(iOp).args = c{2 * iOp};
end
end

% EOF
