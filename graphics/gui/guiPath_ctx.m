function ctx = guiPath_ctx(basepath, basename)
% GUIPATH_CTX Shared context (with a file cache) for guiPath_src / guiPath_load.
%
%   ctx = GUIPATH_CTX(basepath, basename) returns a struct carrying the session
%   location plus a containers.Map cache so a file / computed source is loaded at
%   most once while a whole cfgData is materialized (and reused if the same ctx is
%   passed to later materialize calls, e.g. across preset switches).
%
%   See also guiPath_load, guiPath_src.
%
%   HISTORY:
%       Created: 05 Jul 2026 - declarative redesign.

ctx = struct('basepath', basepath, 'basename', basename, ...
    'cache', containers.Map('KeyType', 'char', 'ValueType', 'any'));
end
