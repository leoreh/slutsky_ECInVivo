function SU = getSU(basepath, spikes, saveVar)

% receives spikes struct and returns same struct but only SU
%
% INPUT:
%   basepath            path to recording folder {pwd}.
%   spikes              struct of spikes (see getSpikes)
%
% CALLS:
%   getSpikes (optional)
%
% 18 dec 18 LH and RA.  Updates:
% 23 dec 18 LH          saveVar

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
if nargs < 3 || isempty(saveVar)
    saveVar = true;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fieldSpikes = fieldnames(spikes);
nclu = length(spikes.UID);
c = cell(length(fieldSpikes), 1);
SU = cell2struct(c, fieldSpikes);
SU = spikes;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Go over fields
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k = 1;
for i = 1 : length(fieldSpikes)
    f = fieldSpikes{i};
    if length(spikes.(f)) <= 1 || isa(spikes.(f), 'char')
        SU.(f) = spikes.(f);
        continue
    end
    for j = 1 : nclu
        if spikes.su(j) == 0
            if isa(spikes.(f), 'double')
                SU.(f)(j) = -999;
            elseif isa(spikes.(f), 'cell')
                SU.(f){j} = -999;
            end
            k = k + 1;
        end
    end
    SU.(f)(SU.(f) == -999) = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save fet
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if saveVar
    [~, filename, ~] = fileparts(basepath);
    save([filename, '.SU'], 'SU')
end

end

% EOF