function units = selectUnits(varargin)

% selects specific units of a session based on multiple params. arranges
% the units in two logical rows (rs; fs) and saves it in the session info
% struct (ce format) and as a separate file. 
%
% INPUT:
%   basepath        char. path to session folder {pwd}
%   spikes          struct
%   cm              struct (cell_metrics)
%   fr              struct. see calc_fr.m 
%   grp             numeric. spike groups to include 
%   frBoundries     2 x 2 mat. include only units with mfr within
%                   frBoundries. 1st row rs; 2nd row fs
%   forceA          logical. re-select units even if units var exists
%   saveVar         logical 
%   altClean        numeric. alternative selection of criteria to define
%                   clean units.
%   graphics        logical
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
addParameter(p, 'forceA', false, @islogical);
addParameter(p, 'saveVar', true, @islogical);
addParameter(p, 'graphics', true, @islogical);
addParameter(p, 'altClean', 1, @isnumeric);

parse(p, varargin{:})
basepath        = p.Results.basepath;
spikes          = p.Results.spikes;
cm              = p.Results.cm;
fr              = p.Results.fr;
grp             = p.Results.grp;
frBoundries     = p.Results.frBoundries;
forceA          = p.Results.forceA;
saveVar         = p.Results.saveVar;
graphics        = p.Results.graphics;
altClean        = p.Results.altClean;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% files
cd(basepath)
[~, basename] = fileparts(basepath);
spktimesfile = fullfile(basepath, [basename, '.spktimes.mat']);
spikesfile = fullfile(basepath, [basename, '.spikes.cellinfo.mat']);
cmfile = fullfile(basepath, [basename, '.cell_metrics.cellinfo.mat']);
frfile = fullfile(basepath, [basename, '.fr.mat']);
unitsfile = fullfile(basepath, [basename, '.units.mat']);

if exist(unitsfile, 'file') && ~forceA
    load(unitsfile, 'units')
    return
end

if isempty(spikes)
    if exist(spikesfile, 'file')
        load(spikesfile)
    end
end
ngrps = length(unique(spikes.shankID));

if isempty(cm)
    load(cmfile, 'cell_metrics')
    cm = cell_metrics;
end
if isempty(fr)
    load(frfile)
end

nunits = length(fr.mfr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% organize
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% su vs mu
units.su = ones(nunits, 1);    % override
if isfield(spikes, 'su')
    units.su = spikes.su';
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
units.grp = grpidx;
% check cell explorer shank id
if any(spikes.shankID ~= sort(spikes.shankID))
    error('Possible error in loadSpikes')
end


% cell class
units.rs = strcmp(cm.putativeCellType, 'Pyramidal Cell');
units.wide = strcmp(cm.putativeCellType, 'Wide Interneuron');
units.fs = strcmp(cm.putativeCellType, 'Narrow Interneuron');

% fr in boundries during baseline
if isempty(frBoundries)
    frBoundries = [0 Inf; 0 Inf];
end
mfrRS = units.rs' & fr.mfr > frBoundries(1, 1) & fr.mfr < frBoundries(1, 2);
mfrFS = units.fs' & fr.mfr > frBoundries(2, 1) & fr.mfr < frBoundries(2, 2);
units.mfrBL = mfrRS | mfrFS;

% fr continuous during baseline and throughout recording
bl_idx = fr.tstamps > fr.info.winBL(1) & fr.tstamps < fr.info.winBL(2);
cnt_bl = sum(fr.strd(:, bl_idx) > frBoundries(1, 1), 2) / sum(bl_idx) > 0.5;
cnt_all = sum(fr.strd > frBoundries(1, 1), 2) / size(fr.strd, 2) > 0.5;
units.cnt = cnt_bl & cnt_all;

% fr "stable" during baseline
units.stable = fr.mfr_std < fr.mfr;
units.stable = ones(nunits, 1);    % override

% fr gini coeff
units.gini = fr.gini_unit <= 0.5 | isnan(fr.gini_unit);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% choose what is clean and save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% select criteria
switch altClean
    case 1 % all
        units.clean(1, :) = logical(units.rs & units.su' & units.grp &...
            units.gini' & units.stable' & units.mfrBL' & units.cnt');
        units.clean(2, :) = logical(units.fs & units.su' & units.grp &...
            units.gini' & units.stable' & units.mfrBL' & units.cnt');

    case 2 % none
        units.clean(1, :) = logical(units.rs) & units.grp;
        units.clean(2, :) = logical(units.fs) & units.grp;

    case 3 % all but stability
        units.clean(1, :) = logical(units.rs & units.su' & units.grp &...
            units.mfrBL');
        units.clean(2, :) = logical(units.fs & units.su' & units.grp &...
            units.mfrBL');
end

% info
units.info.frBoundries = frBoundries;
units.info.grp = grp;

% nunits organized [rs, fs, other (rows); grp1, grp2,... grpn (columns)]
units.nunits = nan(3, ngrps);
for igrp = 1 : ngrps
    grpidx = spikes.shankID == igrp;
    
    for unitType = [1, 2]
        units.nunits(unitType, igrp) = sum(units.clean(unitType, :) & grpidx);
    end
    units.nunits(3, igrp) = sum(~any(units.clean) & grpidx);

end

% save
if saveVar
    save(unitsfile, 'units')
end

% graphics
if graphics
    plot_nunits('basepaths', {basepath}, 'saveFig', true)
end

end

% EOF