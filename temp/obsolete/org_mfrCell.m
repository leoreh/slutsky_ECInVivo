function [mfrCell, gainCell] = org_mfrCell(varargin)

% organizes mfr of units / tetrodes according to states and other params
% (e.g. spiking group, fr boundries, etc). can also organize gain factor
% (nrem - wake).
%
% INPUT: 
%   spikes          struct. CE format.
%   cm              struct. cell metrics from CE. 
%   fr              struct. can be sr / fr. see firingRate.m
%   timebins        numeric mat of n x 2. each row is the start and end of
%                   each time bin [s].
%   stateIdx        numeric vec. indices of states to include in the cell
%                   array. e.g. [1, 4] will include WAKE and NREM. if empty
%                   then will take average fr / sr regardless of state.
%   dataType        string. can be 'mu' and then array will include
%                   spike rate, or 'su and then array will include firing
%                   rate of pyr and int
%   grp             numeric vec. selected tetrodes (spike groups).
%   suFlag          logical. include only well-isolated su or all.
%                   determined by spikes.su
%   frBoundries     numeric 2 x 1. include only units with mean fr in these
%                   boundries
%
% OUTPUT
%   mfrCell         array where each cell contains a mat of mean firing
%                   rates per unit or tetrode (row) in each time bin
%                   (column). if mu than mfrCell is 1d with nstates
%                   columns. if su than mfrCell has two rows, one for
%                   pyr and one for int.
%   gainCell        array where each cell contains a mat of gain fire per
%                   unit (row) in each time bin (column). gain cell has two
%                   rows, one for pyr and one for int. if dataType = 'mu'
%                   or stateIdx is empty, will return empty cell.
%
% DEPENDENCIES
%
% TO DO LIST
%
% 24 oct 21 LH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addOptional(p, 'spikes', []);
addOptional(p, 'cm', []);
addOptional(p, 'fr', []);
addOptional(p, 'dataType', 'mu', @ischar);
addOptional(p, 'timebins', [], @isnumeric);
addOptional(p, 'stateIdx', [], @isnumeric);
addOptional(p, 'grp', [], @isnumeric);
addOptional(p, 'frBoundries', [], @isnumeric);
addOptional(p, 'suFlag', [], @islogical);

parse(p, varargin{:})
spikes = p.Results.spikes;
cm = p.Results.cm;
fr = p.Results.fr;
dataType = p.Results.dataType;
timebins = p.Results.timebins;
stateIdx = p.Results.stateIdx;
grp = p.Results.grp;
frBoundries = p.Results.frBoundries;
suFlag = p.Results.suFlag;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize
gainCell = cell(2, 1);

if strcmp(dataType, 'su')
    if ~isempty(spikes) && ~isempty(cm) && isfield(spikes, 'su')
    pyrUnits = selectUnits(spikes, cm, fr, suFlag, grp, frBoundries, 'pyr');
    intUnits = selectUnits(spikes, cm, fr, suFlag, grp, frBoundries, 'int');
    else
        pyrUnits = ones(1, size(fr.strd, 1));
        intUnits = ones(1, size(fr.strd, 1));
        warning('\n\nSome input may be missing. selecting all units\n\n')
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% organize mfr cell
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isempty(stateIdx)
    for istate = 1 : length(stateIdx)
        for ibin = 1 : size(timebins, 1)
            
            % tstamps relavent to timebin
            tidx = fr.states.tstamps{stateIdx(istate)} > timebins(ibin, 1) &...
                fr.states.tstamps{stateIdx(istate)} < timebins(ibin, 2);
            
            % firing rate data
            switch dataType
                case 'mu'
                    mfrCell{istate}(:, ibin) =...
                        mean(fr.states.fr{stateIdx(istate)}(grp, tidx), 2, 'omitnan');
                    
                case 'su'                    
                    mfrCell{1, istate}(:, ibin) =...         % pyr
                        mean(fr.states.fr{stateIdx(istate)}(pyrUnits, tidx), 2, 'omitnan');
                    mfrCell{2, istate}(:, ibin) =...         % int
                        mean(fr.states.fr{stateIdx(istate)}(intUnits, tidx), 2, 'omitnan');

            end
        end
    end
else
    % regardless of state
    for ibin = 1 : size(timebins, 1)
        
        % tstamps relavent to timebin
        tidx = fr.tstamps > timebins(ibin, 1) &...
            fr.tstamps < timebins(ibin, 2);
        
        % firing rate data
        switch dataType
            case 'mu'
                mfrCell(:, ibin) = mean(fr.strd(grp, tidx), 2, 'omitnan');
                
            case 'su'
                mfrCell{1}(:, ibin) = mean(fr.strd(pyrUnits, tidx), 2, 'omitnan');               
                mfrCell{2}(:, ibin) = mean(fr.strd(intUnits, tidx), 2, 'omitnan');
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% organize gain factor cell
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      

for ibin = 1 : size(timebins, 1)
    
    % tstamps relavent to timebin (wake)
    tidxWake = fr.states.tstamps{1} > timebins(ibin, 1) &...
        fr.states.tstamps{1} < timebins(ibin, 2);
    tidxNREM = fr.states.tstamps{4} > timebins(ibin, 1) &...
        fr.states.tstamps{4} < timebins(ibin, 2);
    
    % rs data
    wakeMat = fr.states.fr{1}(pyrUnits, tidxWake);
    nremMat = fr.states.fr{4}(pyrUnits, tidxNREM);
    gainCell{1}(:, ibin) = (mean(nremMat') - mean(wakeMat')) ./ ...
        max([max(wakeMat'); max(nremMat')]) * 100;
    
    % fs data
    wakeMat = fr.states.fr{1}(intUnits, tidxWake);
    nremMat = fr.states.fr{4}(intUnits, tidxNREM);
    gainCell{2}(:, ibin) = (mean(nremMat') - mean(wakeMat')) ./ ...
        max([max(wakeMat'); max(nremMat')]) * 100;
end
 
end

% EOF