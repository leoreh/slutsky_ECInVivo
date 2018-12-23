function spikeSU = getSU(basepath, spikes)

% receives spikes struct and returns same struct but only SU
%
% INPUT:
%   basepath    path to recording folder {pwd}.
%   spikes      struct of spikes (see getSpikes)
%
% CALLS:
%   getSpikes (optional)
%
% TO DO LIST:
%
% 18 dec 18 LH and RA.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nargs = nargin;
if nargs < 1 || isempty(basepath)
    basepath = pwd;
end
if nargs < 2 || isempty(spikes)
    warning('spikes will be loaded from %s', basepath)
    spikes = getSpikes('basepath', basepath);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fieldSpikes = fieldnames(spikes);
nclu = length(spikes.UID);
c = cell(length(fieldSpikes), 1);
spikesSU = cell2struct(c, fieldSpikes);
spikesSU = spikes;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Go over fields
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k = 1;
for i = 1 : length(fieldSpikes)
    f = fieldSpikes{i};
    if length(spikes.(f)) <= 1 || isa(spikes.(f), 'char')
        spikesSU.(f) = spikes.(f);
        continue
    end
    for j = 1 : nclu
        if spikes.su(j) == 0
            if isa(spikes.(f), 'double')
                spikesSU.(f)(j) = -999;
            elseif isa(spikes.(f), 'cell')
                spikesSU.s(f){j} = -999;
            end
            k = k + 1;
        end
    end
    spikesSU.(f)(spikesSU.(f) == -999) = [];
end

end

% EOF