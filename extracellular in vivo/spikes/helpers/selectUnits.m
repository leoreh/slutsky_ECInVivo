function units = selectUnits(spikes, cm, fr, suFlag, grp, frBoundries, unitClass)
    
    % su vs mu
    su = ones(length(spikes.times), 1);    % override
    if isfield(spikes, 'su') && suFlag
        su = spikes.su';
    end
    
    % tetrode
    grpidx = zeros(1, length(spikes.shankID));
    for ii = 1 : length(grp)
        grpidx = grpidx | spikes.shankID == grp(ii);
    end
       
    % cell class
    pyr = strcmp(cm.putativeCellType, 'Pyramidal Cell');
    wide = strcmp(cm.putativeCellType, 'Wide Interneuron');
    int = strcmp(cm.putativeCellType, 'Narrow Interneuron');
    
    % mfr
    mfrRS = pyr' & fr.mfr > frBoundries(1, 1) & fr.mfr < frBoundries(1, 2);
    mfrFS = int' | wide' & fr.mfr > frBoundries(2, 1) & fr.mfr < frBoundries(2, 2);
    mfrStable = [std(fr.strd') < mean(fr.strd')]';    
    mfrunits = mfrStable & (mfrRS | mfrFS);
    
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
end

% EOF