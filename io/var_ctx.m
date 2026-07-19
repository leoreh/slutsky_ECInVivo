function ctx = var_ctx(basepath, basename)

% Shared context (a file cache) for var_load / var_fetch.
%
% Carries the session location plus a containers.Map so a file or binary is
% read at most once while a varMap is filled, and reused if the same ctx is
% passed to later calls (e.g. across preset switches in a viewer).
%
% INPUTS
% - basepath        <char> session folder.
% - basename        <char>(opt) file stem. Default: folder name of basepath.
%
% OUTPUTS
% - ctx             <struct> .basepath .basename .cache (containers.Map).
%
% HISTORY
% - 260719          created (unified var_* I/O layer; was guiPath_ctx).

if nargin < 1 || isempty(basepath), basepath = pwd; end
if nargin < 2 || isempty(basename), [~, basename] = fileparts(basepath); end

ctx = struct('basepath', basepath, 'basename', basename, ...
    'cache', containers.Map('KeyType', 'char', 'ValueType', 'any'));

end

% EOF
