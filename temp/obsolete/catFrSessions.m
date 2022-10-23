function frCat = catFrSessions(varargin)

% reveives a cell of basepaths, loads the fr struct, and combines fields of
% interest (e.g. state ratio) across units. also returns a logical mat of
% units based on criteria. can load a specific time bin from frBins instead
% of fr.
% 
% INPUT
%   basepaths   cell of char to recording sessions
%   binidx      numeric. idx to one of the time bins in frBins. if exist
%               then frBins(binidx) will be used instead of fr
% 
% OUTPUT
%   frCat       struct
%      
% 23 jan 22 LH      

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;

p = inputParser;
addOptional(p, 'basepaths', @iscell);
addOptional(p, 'binidx', [], @isnumeric);

parse(p, varargin{:})
basepaths   = p.Results.basepaths;
binidx      = p.Results.binidx;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% params
stateIdx = [1 : 5];
frBoundries = [0 Inf; 0 Inf];
stableFlag = false;
giniFlag = false;
  
% load vars from each session
varsFile = ["fr"; "units"; "fr_bins"];
varsName = ["fr"; "units"; "frBins"];
v = getSessionVars('basepaths', basepaths, 'varsFile', varsFile,...
    'varsName', varsName);
nsessions = length(basepaths);

% if binidx was selected, replace the fr struct with one of the timebins
% in frBins
if binidx
    for isession = 1 : length(v)
        v(isession).fr = v(isession).frBins(binidx);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize
stateMfr = [];
stateRat = [];
units = [];
unitsGini = [];
unitsMfr = [];
mfr = [];
giniCoeff = [];
stableIdx = [];
sessionIdx = [];
tstamps = [];

for isession = 1 : nsessions
    basepath = basepaths{isession};
    cd(basepath)    
    
    % mean firing rate
    mfr = [mfr; v(isession).fr.mfr];
    
    % tstamps
    tstamps = [tstamps, v(isession).fr.tstamps];
    
    % mean firing rate for each unit in states
    stateMfr_tmp = cellfun(@(x) mean(x, 2), v(isession).fr.states.fr, 'uni', false);
    stateMfr = [stateMfr, cell2mat(stateMfr_tmp(stateIdx))'];
    
    % state ratio (NREM - WAKE)
    stateRat = [stateRat; squeeze(v(isession).fr.states.ratio(4, 1, :))];
    
    % select specific units
    units = logical([units, v(isession).units.idx]);
    unitsGini = logical([unitsGini; v(isession).units.gini]);
    unitsMfr = logical([unitsMfr; v(isession).units.mfrBL]);
    nunits = size(v(isession).units.idx, 2);
    sessionIdx = [sessionIdx; ones(nunits, 1) * isession];
    
    % gini and stable
%     giniCoeff = [giniCoeff; v(isession).fr.gini_unit];
    stableIdx = logical([stableIdx; v(isession).fr.stable]);
end

% organize in struct for output
frCat.stateMfr = stateMfr;
frCat.stateRat = stateRat;
frCat.units = units;
frCat.unitsGini = unitsGini;
frCat.unitsMfr = unitsMfr;
frCat.mfr = mfr;
% frCat.giniCoeff = giniCoeff;
frCat.stableIdx = stableIdx;
frCat.sessionIdx = sessionIdx;
frCat.tstamps = tstamps;

end

% EOF