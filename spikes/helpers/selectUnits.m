function units = selectUnits(varargin)

% selects specific units of a session based on multiple params. arranges
% the units in two logical rows (rs; fs) and saves it in the session info
% struct (ce format) and as a separate file. 
%
% INPUT:
%   basepath        char. path to session folder {pwd}
%   spikes          struct
%   cm              struct (cell_metrics)
%   fr              struct. see firingRate.m 
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
addParameter(p, 'forceA', false, @islogical);
addParameter(p, 'saveVar', true, @islogical);

parse(p, varargin{:})
basepath        = p.Results.basepath;
spikes          = p.Results.spikes;
cm              = p.Results.cm;
fr              = p.Results.fr;
grp             = p.Results.grp;
frBoundries     = p.Results.frBoundries;
forceA          = p.Results.forceA;
saveVar         = p.Results.saveVar;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% files
cd(basepath)
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
    if exist(spikesfile, 'file')
        load(spikesfile)
    end
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
if isfield(spikes, 'su')
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
rs = strcmp(cm.putativeCellType, 'Pyramidal Cell');
wide = strcmp(cm.putativeCellType, 'Wide Interneuron');
fs = strcmp(cm.putativeCellType, 'Narrow Interneuron');

% fr in boundries during baseline
if isempty(frBoundries)
    frBoundries = [0 Inf; 0 Inf];
end
mfrRS = rs' & fr.mfr > frBoundries(1, 1) & fr.mfr < frBoundries(1, 2);
mfrFS = fs' & fr.mfr > frBoundries(2, 1) & fr.mfr < frBoundries(2, 2);

% fr continuous during baseline and throughout recording
bl_idx = fr.tstamps > fr.info.winBL(1) & fr.tstamps < fr.info.winBL(2);
cnt_bl = sum(fr.strd(:, bl_idx) > 0.01, 2) / sum(bl_idx) > 0.8;
cnt_all = sum(fr.strd > 0.01, 2) / size(fr.strd, 2) > 0.8;

% fr "stable" during baseline
fr.stable = fr.mfr_std < fr.mfr;

% combine
unitsidx(1, :) = rs & su' & grpidx;
unitsidx(2, :) = fs & su' & grpidx;

% arrange struct
units.idx = logical(unitsidx);
units.mfrBL = mfrRS | mfrFS;
units.cnt = cnt_bl & cnt_all;
units.stable = fr.stable;
units.gini = fr.gini_unit <= 0.5 | isnan(fr.gini_unit);
units.grp = grpidx;
units.rs = rs;
units.fs = fs;
units.wide = wide;
units.su = su;
units.info.frBoundries = frBoundries;
units.info.grp = grp;

% save
if saveVar
    save(unitsfile, 'units')
end

end

% EOF