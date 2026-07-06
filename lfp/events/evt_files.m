function files = evt_files(basepath, basename, name)
% EVT_FILES Standard output file paths for an event modality.
%
%   files = EVT_FILES(basepath, basename, name)
%
%   SUMMARY:
%       Shared filename builder for both event pipelines. Returns the canonical
%       output paths for a modality, keyed by its short name ('ripp' or 'ed'),
%       so the wrappers and their readers agree on one naming scheme. The phase
%       file is ripples-only and is added by that wrapper; the per-bout states
%       file is produced by evt_states itself and is not returned here.
%
%   INPUTS:
%       basepath - (Char) Session directory.
%       basename - (Char) File stem (may differ from the folder name).
%       name     - (Char) Modality tag, 'ripp' or 'ed'.
%
%   OUTPUTS:
%       files    - (Struct) Absolute paths:
%           .evt     - <basename>.<name>.mat      (events, curation)
%           .maps    - <basename>.<name>Maps.mat  (per-event LFP maps)
%           .spks    - <basename>.<name>Spks.mat  (per-unit stats + PETH, light)
%           .spkMaps - <basename>.<name>SpkMaps.mat (3D spike raster, heavy)
%
%   DEPENDENCIES:
%       None.
%
%   HISTORY:
%       Created: 260706 (unifies the wrappers' output-path construction).

files.evt     = fullfile(basepath, [basename, '.', name, '.mat']);
files.maps    = fullfile(basepath, [basename, '.', name, 'Maps.mat']);
files.spks    = fullfile(basepath, [basename, '.', name, 'Spks.mat']);
files.spkMaps = fullfile(basepath, [basename, '.', name, 'SpkMaps.mat']);

end     % EOF
