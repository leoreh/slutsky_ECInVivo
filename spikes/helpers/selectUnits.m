function units = selectUnits(spikes, cm, fr, suFlag, grp, frBoundries, unitClass)

% selects specific units of a session based on multiple params

if isempty(frBoundries)
    frBoundries = [0 Inf; 0 Inf];
end


nunits = length(fr.mfr);

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

% mfr
mfrRS = pyr' & fr.mfr > frBoundries(1, 1) & fr.mfr < frBoundries(1, 2);
mfrFS = int' | wide' & fr.mfr > frBoundries(2, 1) & fr.mfr < frBoundries(2, 2);
mfrunits = fr.stable & fr.bl_thr & (mfrRS | mfrFS);

% combine
if strcmp(unitClass, 'pyr')
    units = pyr & su' & grpidx & mfrunits';
elseif strcmp(unitClass, 'int')
    units = int & su' & grpidx & mfrunits';
elseif strcmp(unitClass, 'wide')
    units = wide & su' & grpidx & mfrunits';
else
    units = su' & grpidx & mfrunits';     % override
end
units = logical(units(:));

end

% EOF