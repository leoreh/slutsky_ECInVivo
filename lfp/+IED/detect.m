function ied = detect(varargin)
% detects inter-ictal discharges from LFP.
%
% INPUT (in 1st position):
%   IED.data object, which already contain all the name values
% OR:
% INPUT (name-value):
%   sig         signal for detection
%   fs          sampling frequency
%   thr         vector of two elements. first is thr in [z-scores] and
%               second is thr in [mV]. if one specified than the other is
%               merely calculated. if both specified than both used for
%               detection. values should be positive even if thrDir is
%               negative
%   thrDir      string. direction of threshold detection. can be
%               {'positive'}, 'negative', or 'both
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
switch ied.thrDir
    case 'positive'
        if ied.thr(1) == 0
            thresholded = ied.sig > ied.thr(2);
            % find thr in z scores
            ied.thr(1) = (ied.thr(2) - mean(ied.sig)) / std(ied.sig);
        elseif ied.thr(2) == 0
            thresholded = zscore(ied.sig) > ied.thr(1);
            % find thr in mV
            ied.thr(2) = mean(ied.sig) + ied.thr(1) * std(ied.sig);
        else
            thresholded = zscore(ied.sig) > ied.thr(1) & ied.sig > ied.thr(2);
        end
    case 'negative'
        if ied.thr(1) == 0
            thresholded = ied.sig < -ied.thr(2);
            % find thr in z scores
            ied.thr(1) = -(ied.thr(2) - mean(ied.sig)) / std(ied.sig);
        elseif ied.thr(2) == 0
            thresholded = zscore(ied.sig) < -ied.thr(1);
            % find thr in mV
            ied.thr(2) = -(mean(ied.sig) + ied.thr(1) * std(ied.sig));
        else
            thresholded = zscore(ied.sig) < -ied.thr(1) & ied.sig < -ied.thr(2);
        end
    case 'both'
        if ied.thr(1) == 0
            thresholded = ied.sig > ied.thr(2) | ied.sig < -ied.thr(2);
            % find thr in z scores
            ied.thr(1) = (ied.thr(2) - mean(ied.sig)) / std(ied.sig);
        elseif ied.thr(2) == 0
            thresholded = zscore(ied.sig) > ied.thr(1) | zscore(ied.sig) < -ied.thr(1);
            % find thr in mV
            ied.thr(2) = mean(ied.sig) + ied.thr(1) * std(ied.sig);
        else
            thresholded = (zscore(ied.sig) > ied.thr(1) & ied.sig > ied.thr(2)) |...
                (zscore(ied.sig) < -ied.thr(1) & ied.sig < -ied.thr(2));
        end
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
[pos,idx] = uniquetol(pos,interDur,'lowest','DataScale',1);
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
% baseline = mean(ied.sig) + 1 * std(ied.sig);
for i = 1 : nspks
    
    % check that peak-to-peak in narrow window is greater than thr in mV
    win = ied.sig(pos(i) - round(0.015 * ied.fs) : pos(i) + round(0.015 * ied.fs));
    amp(i) = peak2peak(win);
    if  amp(i) < ied.thr(2)
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