function drft = drift_calc(fr_mat, tstamps, varargin)

% core function to calculate drift. recieves a matrix of firing rates
% across time, and a vector of timestamps. divides the recording to
% windows, calculating a population vector (PV) per window, the pairwise
% correlation between PVs, and the drift rate via linear fit.
%
% INPUT:
%   fr_mat          numeric mat of firing rates [units x tstamps]
%   tstamps         numeric vec of timestamps, relative to recording onset
%   winsize         numeric. window size for calculating PVs [s] {3600}
%   nwin            numeric. no of windows to divide the recording. will
%                   override winsize
%   thrLin          numeric. when calc lin fit, include only win with thrLin corr pairs 
%   thrFr           numeric. only units with fr > thrRr are included in the PV
%   thrWin          numeric. only windows calculated with more bins then thrWin are included
%   limUnit         numeric. randomely select limUnit units in each PV 
%   basepath        string. path to recording folder {pwd}
%   graphics        logical. plot {false}
%
% DEPENDENCIES:
%   drift_plot
%
% TO DO LIST:
%   * use spktimes instead of pre-calculated firing rates. this is especially
%   important for the fr criterion. perhaps it would be best to calculate
%   FRs in 1 hour bins before this function (e.g., in drift_file)
%   * make sure there is enough data points per window (done)
%
% 22 may 24 LH      based on Lee's code

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addOptional(p, 'basepath', pwd);
addOptional(p, 'winsize', 3600, @isnumeric);
addOptional(p, 'nwin', [], @isnumeric);
addOptional(p, 'thrLin', [], @isnumeric);
addOptional(p, 'thrFr', [], @isnumeric);
addOptional(p, 'thrWin', [], @isnumeric);
addOptional(p, 'limUnit', [], @isnumeric);
addOptional(p, 'graphics', false, @islogical);

parse(p, varargin{:})
basepath        = p.Results.basepath;
winsize         = p.Results.winsize;
nwin            = p.Results.nwin;
thrLin          = p.Results.thrLin;
thrFr           = p.Results.thrFr;
thrWin          = p.Results.thrWin;
limUnit         = p.Results.limUnit;
graphics        = p.Results.graphics;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% file params
[~, basename] = fileparts(basepath);

% recording params
recLen = tstamps(end);
nunits = size(fr_mat, 1);

% divide the recording to nwin (alt 1) or to wins of precisley 1 hr (alt 2)
if isempty(nwin)
    wins = n2chunks('n', recLen, 'chunksize', winsize);
    wins(end, :) = [];
else
    wins = n2chunks('n', recLen, 'nchunks', nwin);
end
nwin = size(wins, 1);
ncorr = nwin - 1;

% prepare output
drft = struct('dt_corr', [], 'm_corr', [], 'lin_coef', [], 'drate', []);
drft.info = struct('winsize', [], 'thrLin', [], 'thrIdx', [], 'thrFr', [],...
'thrWin', [], 'rmUnits', [], 'rmWins', []);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% analyze
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calc PVs - population fr vectors for each window
pv = nan(nunits, nwin);
for iwin = 1 : nwin
    winIdx = InIntervals(tstamps, wins(iwin, :));
    npnts(iwin) = sum(winIdx);
    fr_vec = mean(fr_mat(:, winIdx), 2);
    pv(:, iwin) = fr_vec / norm(fr_vec);
end

% apply criterion such that only units with fr > thrFr are included in the
% PV
rmUnits = min(pv, [], 2) < thrFr;
pv(rmUnits, :) = [];
nunits = size(pv, 1);

% randomely select a subset of units
if ~isempty(limUnit)
    if limUnit <= nunits
        unitIdx = randsample([1 : nunits], limUnit);
        pv = pv(unitIdx, :);
    else
        warning('limUnit smaller then nunits')
        return
    end
end

% apply criterion such that only windows calculated with more bins then
% thrWin are included
rmWins = npnts < thrWin;
pv(:, rmWins) = nan(nunits, sum(rmWins));

% calc correlation between pv pairs, and organize by delta time
pv_corr = corr(pv);
dt_corr = cell(ncorr, 1);
for iwin = 1 : ncorr
    dt_corr{iwin} = diag(pv_corr, iwin);
end
dt_corr = cell2padmat(dt_corr, 2)';
m_corr = mean(dt_corr, 2, 'omitnan');

% limit wins used for linear fit to those with more than linThr
% correlations pairs
if ~isempty(thrLin)
    thrIdx = sum(~isnan(dt_corr)) >= thrLin;
else
    thrIdx = true(1, ncorr);
end

% exclude DTs that are all nan due to any of the previous criterion
thrIdx = thrIdx & ~isnan(m_corr)';

% calc linear fit to pv correlations, where slope is drift rate
xaxis = [1 : ncorr];
lin_coef = polyfit(xaxis(thrIdx), m_corr(thrIdx), 1);
drate = abs(lin_coef(1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% organize output and plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

drft.dt_corr = dt_corr;
drft.m_corr = m_corr;
drft.lin_coef = lin_coef;
drft.drate = drate;
drft.info.winsize = winsize;
drft.info.thrLin = thrLin;
drft.info.thrIdx = thrIdx;
drft.info.thrFr = thrFr;
drft.info.thrWin = thrWin;
drft.info.rmUnits = rmUnits;
drft.info.rmWins = rmWins;


if graphics
    drift_plot(drft)
end

end

% EOF

