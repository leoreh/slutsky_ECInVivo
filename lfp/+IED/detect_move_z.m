function ied = detect_move_z(varargin)
% detects inter-ictal discharges from LFP, using a moving z-score threshold
%
% INPUT (in 1st position):
%   IED.data object, which already contain all the name values
% OR:
% INPUT (name-value, required):
%   sig         signal for detection
%   fs          sampling frequency
%   thr         positive scalar, z-score threshold to pass for moving-z
%               threshold
%
% INPUT (name-value, optional)
%   thrDir      text scalar. direction of threshold detection. can be
%               'positive', 'negative', or {'both'}
%   emg         emg signal. Note that it MUST have the same sampling rate
%               as sig, and be of the same length.
%
% OUTPUT
%   IED.data object
%
% Based on getIIS by LH (see +IED/legacy folder)
% By: LdM 
% Published: 230827
%
%   see also IED.data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check if input is already IED.data object
if ~isa(varargin{1},'IED.data')
    % place everything into ied
    ied = IED.data(varargin{:});
else
    % use given object
    ied = varargin{1};
end

% params
interDur = round(ied.fs * 0.025);       % samples
lowthr = 0.2;                       % mV

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% detect discharges
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% threshold
calc_win = ied.fs*5; % 5 sec
sig_mu = movmean(ied.sig,calc_win);
sig_sigma = movstd(ied.sig,calc_win);
local_z = (ied.sig - sig_mu)./sig_sigma;
switch ied.thrDir
    case 'positive'
        thresholded = local_z > ied.thr;
    case 'negative'
        thresholded = local_z < ied.thr;
    case 'both'
        thresholded = (local_z > ied.thr) | (local_z < -ied.thr);
end
thresholded = thresholded(:); % force col vec
iie = find([0; diff(thresholded) > 0]);

% return if no IIS are detected
if isempty(iie)
    fprintf('\nno inter-ictal discharges detected\n');
    return
end

% remove crossings that are too close together. this effectively limits the
% maximum IIS burst frequency to 1 / interDur. the division of interDur is
% so that the search for max later on extends beyond the deleted
% crossings
ii = find(diff(iie) < interDur / 2);
while ~isempty(ii)
    iie(ii + 1) = [];
    ii = find(diff(iie) < interDur / 2);
end

% select local maximum and clip discharges 
peak = zeros(length(iie), 1);
pos = zeros(length(iie), 1);
rmidx = [];
for i = 1 : length(iie)
    % remove last / first iis if incomplete
    if iie(i) + interDur > length(ied.sig) || iie(i) - interDur < 1
        rmidx = [rmidx, i];
        continue
    end
    % adjust crossing to local max / min and then clip
    localseg = ied.sig(iie(i) - interDur : iie(i) + interDur);
    if abs(max(localseg)) > abs(min(localseg))
        [peak(i), pos(i)] = max(localseg);
    else
        [peak(i), pos(i)] = min(localseg);
    end
    pos(i) = iie(i) - interDur + pos(i) -1; %-1 as the "+" make you count 1 sample twice
end
peak(rmidx) = [];
pos(rmidx) = [];

% due to overlapping localseg, 2 positions may be equal or out of order.
% fix that.
[pos,idx] = uniquetol(pos,interDur,'DataScale',1);
% [pos,idx] = unique(pos);
peak = peak(idx);

% if there is no trough between peaks select highest peak. this is
% instead of demanding for a refractory period via interDur
rmidx = [];
for i = 1 : length(pos) - 1
    low = min(ied.sig(pos(i) : pos(i + 1)));
    if low > lowthr
        rmidx = [rmidx; i + (peak(i) > peak(i + 1))];
    end
end
pos(rmidx) = [];
nspks = length(pos);
fprintf('\nafter detection: %d IED\n', nspks);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% template matching
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rmidx = [];
thr_MV =  ied.thr.*sig_sigma + sig_mu;
% baseline = mean(ied.sig) + 1 * std(ied.sig);
for i = 1 : nspks
    
    % check that peak-to-peak in narrow window is greater than thr in mV
    win = ied.sig(pos(i) - round(0.015 * ied.fs) : pos(i) + round(0.015 * ied.fs));
    amp(i) = peak2peak(win);
    if  amp(i) < thr_MV(pos(i))
        rmidx = [rmidx, i];
    end
    
    % check that signal returns to baseline before and after peak
    % half1 = ied.sig(pos(i) : -1 : pos(i) - round(0.02 * fs));
    % half2 = ied.sig(pos(i) : pos(i) + round(0.02 * fs));
    % if min(half1) > baseline || min(half2) > baseline
    %     rmidx = [rmidx, i];
    % end
end
pos(rmidx) = [];
nspks = length(pos);
fprintf('after template matching: %d IED\n', nspks);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% organize output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% place index inside
ied.pos = pos;

% change analysis stage
ied.status = "pre_curation";

% make sure there is accepted status for every discharge
ied.accepted = true(1,numel(ied.pos));

% init last discharge curated
ied.last_mark = 1;

% EOF