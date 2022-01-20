function units = selectUnits(varargin)

% selects specific units of a session based on multiple params. arranges
% the units in two logical rows (rs; fs) and saves it in the session info
% struct (ce format). can also update the spikes struct with the params of
% selection
%
% INPUT:
%   basepath        char. path to session folder {pwd}
%   spikes          struct
%   cm              struct (cell_metrics)
%   fr              struct. see firingRate.m 
%   suFlag          logical. include only well isolated single units 
%   stableFlag      logical. include only units with stable baseline firing
%                   rate (according to fr struct)
%   giniFlag        logical. include only units with gini coeff > 0.5
%   blFlag          logical. include only units that passed the baseline
%                   mfr threshold
%   grp             numeric. spike groups to include 
%   frBoundries     2 x 2 mat. include only units with mfr within
%                   frBoundries. 1st row rs; 2nd row fs
%   forceA          logical. re-select units even if units var exists
%   saveVar         logical 
%
% DEPENDENCIES:
% 
% TO DO LIST:
%
% 06 feb 21 LH  updates:
% 12 jan 22 LH  save in spikes struct

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addParameter(p, 'basepath', pwd, @ischar);
addParameter(p, 'spikes', []);
addParameter(p, 'cm', []);
addParameter(p, 'fr', []);
addParameter(p, 'grp', []);
addParameter(p, 'frBoundries', [0 Inf; 0 Inf], @isnumeric);
addParameter(p, 'suFlag', true, @islogical);
addParameter(p, 'stableFlag', false, @islogical);
addParameter(p, 'giniFlag', false, @islogical);
addParameter(p, 'blFlag', false, @islogical);
addParameter(p, 'forceA', false, @islogical);
addParameter(p, 'saveVar', true, @islogical);

parse(p, varargin{:})
basepath        = p.Results.basepath;
spikes          = p.Results.spikes;
cm              = p.Results.cm;
fr              = p.Results.fr;
grp             = p.Results.grp;
frBoundries     = p.Results.frBoundries;
suFlag          = p.Results.suFlag;
stableFlag      = p.Results.stableFlag;
giniFlag        = p.Results.giniFlag;
blFlag          = p.Results.blFlag;
forceA          = p.Results.forceA;
saveVar         = p.Results.saveVar;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% files
[~, basename] = fileparts(basepath);
spikesfile = fullfile(basepath, [basename, '.spikes.cellinfo.mat']);
cmfile = fullfile(basepath, [basename, '.cell_metrics.cellinfo.mat']);
frfile = fullfile(basepath, [basename, '.fr.mat']);
unitsfile = fullfile(basepath, [basename, '.units.mat']);

if exist(unitsfile, 'file') && ~forceA
    load(unitsfile, 'units')
    return
end

if isempty(spikes)
    load(spikesfile)
end
if isempty(cm)
    load(cmfile, 'cell_metrics')
    cm = cell_metrics;
end
if isempty(fr)
    load(frfile)
end

nunits = length(fr.mfr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% selection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% su vs mu
su = ones(nunits, 1);    % override
if isfield(spikes, 'su') && suFlag
    su = spikes.su';
end

% tetrode
if isempty(grp)
    grpidx = ones(1, nunits);
else
    grpidx = zeros(1, nunits);
    for igrp = 1 : length(grp)
        grpidx = grpidx | spikes.shankID == grp(igrp);
    end
end

% cell class
pyr = strcmp(cm.putativeCellType, 'Pyramidal Cell');
wide = strcmp(cm.putativeCellType, 'Wide Interneuron');
int = strcmp(cm.putativeCellType, 'Narrow Interneuron');

% fr 
if isempty(frBoundries)
    frBoundries = [0 Inf; 0 Inf];
end
mfrRS = pyr' & fr.mfr > frBoundries(1, 1) & fr.mfr < frBoundries(1, 2);
mfrFS = (int' | wide') & fr.mfr > frBoundries(2, 1) & fr.mfr < frBoundries(2, 2);

mfrStable = ones(nunits, 1);
if stableFlag
    mfrStable = fr.stable;
end

mfrGini = ones(nunits, 1);
if giniFlag
    mfrGini = fr.gini_unit <= 0.5;
end

mfrBL = ones(1, nunits);
if blFlag
    mfrBL = fr.bl_thr;
end

% combine
units(1, :) = pyr & su' & grpidx & mfrRS' & mfrStable' & mfrBL & mfrGini';
units(2, :) = int & su' & grpidx & mfrFS' & mfrStable' & mfrBL & mfrGini';
units = logical(units);

% save
if saveVar
    spikes.units.idx = units;
    spikes.units.frBoundries = frBoundries;
    spikes.units.grp = grp;
    spikes.units.suFlag = suFlag;
    spikes.units.stableFlag = stableFlag;
    spikes.units.blFlag = blFlag;
    save(spikesfile, 'spikes')
    save(unitsfile, 'units')
end

end

% EOF