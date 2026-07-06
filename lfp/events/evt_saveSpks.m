function evt_saveSpks(spks, basepath, basename, name)
% EVT_SAVESPKS Split the spike results into a light stats file + heavy raster.
%
%   EVT_SAVESPKS(spks, basepath, basename, name)
%
%   SUMMARY:
%       Shared spike-output writer for both event pipelines. Splits the evt_spks
%       result by weight so the file the manuscript loads hot stays small:
%           <basename>.<name>SpkMaps.mat - the 3D per-event raster (.su/.mu,
%               each .evt/.ctrl) + .tstamps. Heavy; the source for any PETH.
%           <basename>.<name>Spks.mat    - per-unit stats + per-unit PETH, the
%               same struct without .maps. Light.
%       Each file holds a single variable named for the modality (rippSpkMaps /
%       edSpkMaps and rippSpks / edSpks) so basepaths2vars and load() resolve it
%       by the filename token. Call only when saving; plotting keeps using the
%       full in-memory struct (with .maps), so this writer returns nothing.
%
%   INPUTS:
%       spks     - (Struct) evt_spks result: .maps, .tstamps, per-unit stats,
%                           .peth (and no .events - the wrapper moved it out).
%       basepath - (Char)   Session directory.
%       basename - (Char)   File stem.
%       name     - (Char)   Modality tag, 'ripp' or 'ed'.
%
%   OUTPUTS:
%       None (writes two .mat files).
%
%   DEPENDENCIES:
%       evt_files.
%
%   HISTORY:
%       Created: 260706 (absorbs the duplicated wrapper spike-save split).

files = evt_files(basepath, basename, name);

% Heavy: the 3D raster + its time base, saved as <name>SpkMaps
mapsVar = [name, 'SpkMaps'];
spkMaps = spks.maps;
spkMaps.tstamps = spks.tstamps;
S.(mapsVar) = spkMaps;
save(files.spkMaps, '-struct', 'S', '-v7.3');

% Light: everything except the raster, saved as <name>Spks
lightVar = [name, 'Spks'];
L.(lightVar) = rmfield(spks, 'maps');
save(files.spks, '-struct', 'L', '-v7.3');

end     % EOF
