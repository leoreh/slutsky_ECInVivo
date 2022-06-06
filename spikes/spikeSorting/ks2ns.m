function ks2ns(rez)

% converts KiloSort output (.rez structure) to neurosuite files (fet, res,
% clu, and spk). based in part on Kilosort2Neurosuite (Peterson). extracts
% and detrends waveforms from the dat file (see getSPKfromDat for
% more information). IMPROTANT: XML FILE MUST BE UPDATED TO THE PARAMS USED
% HERE (E.G. SAMPLES PER WAVEFORM)
%
% DEPENDENCIES
%   class2bytes
%   spktimes2ns
%
% TO DO LIST:
%   # find trend based on larger segment without spike
%   # fix clipped spike at end / beginning of recording
%
% 15 jun 20 LH      updates:
% 07 aug 20         fix case where no spikes in grp
% 01 jun 22         ks3

% session params 
session = CE_sessionTemplate(pwd, 'viaGUI', false,...
    'forceDef', true, 'forceL', false, 'saveVar', false);      
basepath = session.general.basePath;
[~, basename] = fileparts(basepath);
nchans = session.extracellular.nChannels;
fs = session.extracellular.sr;
spkgrp = session.extracellular.spikeGroups.channels;
ngrps = length(spkgrp);
nchans_spk = length([spkgrp{:}]);

% organize rez output
st3 = sortrows(rez.st3, [2, 1]);
tmptimes = st3(:, 1);
tmpid = st3(:, 2);
ntmp = length(unique(tmpid));

% find to which spkgrp each template belongs
tmplts = zeros(ntmp, rez.ops.nt0, nchans_spk);
for itmp = 1 : ntmp
    tmplts(itmp, :, :) = [squeeze(rez.U(:, itmp, :)) * squeeze(rez.W(:, itmp, :))']';
    [~, ampMaxCh(itmp)] = max(range(tmplts(itmp, :, :)));
end

% separate spikes according to groups
tmpch = rez.ops.kcoords(ampMaxCh);
for igrp = 1 : ngrps
    tmpgrp = find(tmpch == igrp);
    spkid{igrp} = find(ismember(tmpid, tmpgrp));
    cluid{igrp} = tmpid(spkid{igrp});
    spktimes{igrp} = tmptimes(spkid{igrp});
end

% save spktimes
spkfile = fullfile(basepath, [basename, '.spktimes.mat']);
save(spkfile, 'spktimes')

% create ns files 
for igrp = 1 : ngrps    
    saveNS(cluid{igrp}, 'datatype', 'clu', 'grpid', igrp)
end
spktimes2ns('basepath', basepath, 'fs', fs,...
    'nchans', nchans, 'spkgrp', spkgrp, 'mkClu', false,...
    'dur', [], 't', [], 'grps', [1 : length(spkgrp)],...
    'spkFile', 'temp_wh');


end

% EOF