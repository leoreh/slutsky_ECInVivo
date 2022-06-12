function spikes = cluVal(varargin)

% wrapper for inspecting and validating clusters
% computes:
%   (1) cluster isolation ditance via:
%       (a) L ratio (Schmitzer-Torbert and Redish, 2004)
%       (b) Isolation distance (Harris et al., 2001)
%   (2) ISI contamination
%
% INPUT:
%   basepath        path to recording
%   spikes          struct (see getSpikes)
%   npca            #PCA features for mahal dist {12}.
%                   If sorted with neurosuite should be #spk group * 3
%   spkgrp          array where each cell is the electrodes for a spike group. 
%   mu              vector of predetermined multi units. spikes from these
%                   units will be discarded from distance calculations.
%                   numbering should correspond to spikes.UID
%   force           logical. force recalculation of distance.
%   graphics        logical. plot graphics {true} or not (false)
%   vis             string. show figure {'on'} or not {'off'}
%   saveFig         logical. save figure {1} or not (0)
%   saveVar         logical. save variables (update spikes and save su)
%   session         struct. session info (see CE_sessionTemplate)
%
% OUTPUT:           additional fields to struct (all 1 x nunits):
%   su              SU (1) or MU (0)
%   L
%   iDist
%   ISIratio
%
% DEPENDENCIES:
%   getSpikes       optional (if not given as input)
%   loadNS          load fet
%   cluDist         L ratio and iDist
%   plotCluster     plot waveform and ISI histogram
%
% TO DO LIST:
%   # adapt to cell explorer (done)
%   # recalculate pca if fet files don't exist
%   # adapt UMS2K for cleaning clusters
%
% 03 dec 18 LH      UPDATS: 
% 15 nov 19 LH      separated ISI to standalone function
%                   added predetermined mu
%                   adjusted i\o params
% 25 jun 20 LH      spkgrp and adaptation to cell explorer format
% 04 jul 21 LH      replaced getFet w/ loadNS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addOptional(p, 'spikes', []);
addOptional(p, 'basepath', pwd, @ischar);
addOptional(p, 'npca', 3, @isscalar);
addOptional(p, 'spkgrp', {}, @iscell);
addOptional(p, 'mu', [], @isnumeric);
addOptional(p, 'force', false, @islogical);
addOptional(p, 'graphics', true, @islogical);
addOptional(p, 'vis', 'on', @ischar);
addOptional(p, 'saveFig', true, @islogical);
addOptional(p, 'saveVar', true, @islogical);
addOptional(p, 'session', []);

parse(p, varargin{:})
spikes      = p.Results.spikes;
basepath    = p.Results.basepath;
npca        = p.Results.npca;
spkgrp      = p.Results.spkgrp;
mu          = p.Results.mu;
force       = p.Results.force;
graphics    = p.Results.graphics;
vis         = p.Results.vis;
saveFig     = p.Results.saveFig;
saveVar     = p.Results.saveVar;
session     = p.Results.session;

% sessionInfo
[~, basename] = fileparts(basepath);
if isempty(session)
    sessionName = [basename, '.session.mat'];
    if ~exist(sessionName, 'file')
        session = CE_sessionTemplate(pwd, 'viaGUI', false,...
            'force', true, 'saveVar', true);
    else
        load(sessionName)
    end
end
if isempty(spkgrp)
    spkgrp = session.extracellular.spikeGroups.channels;
end

% load spikes if empty
if isempty(spikes)
    [~, filename] = fileparts(basepath);
    spkname = [filename '.spikes.cellinfo.mat'];
    if exist(spkname, 'file')
        load(spkname)
    else
        error('%s not found', spkname)
    end
end
% make sure spikes has required fields
if ~all(isfield(spikes, {'shankID', 'cluID', 'times'}))
    error('spikes missing required fields')
end

% params
grps = unique(spikes.shankID);
ngrp = length(grps);  % only tetrodes with units
nunits = length(spikes.times);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% waveform params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define a good waveform such that at least two channels exhibit a spike
% amplitude (range) greater than the maximum std on that channel

% count = 1;
% for igrp = 1 : ngrp
%     clu = loadNS('datatype', 'clu', 'session', session, 'grpid', igrp);
%     uclu = unique(clu);
%     spk = loadNS('datatype', 'spk', 'session', session, 'grpid', igrp,...
%         'nspks', length(clu));
%     
%     for iclu = 1 : length(uclu)
%         if uclu(iclu) == 0 || uclu(iclu) == 1
%             continue
%         end
%         wvgrp = spk(:, :, clu == uclu(iclu));
%         
%         goodwv(count) = sum(range([mean(wvgrp, 3)]') > max([std(wvgrp, [], 3)]'));
%         count = count + 1;
%     end
% end    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ISI contamination
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SU = 1% of ISIs < 0.002
ref = 0.002;
thr = 1;
spikes.isiViolation = zeros(1, nunits);
for iunit = 1 : nunits
    nspks(iunit, 1) = length(spikes.times{iunit});
    spikes.isiViolation(1, iunit) = sum(diff(spikes.times{iunit}) < ref) / nspks(iunit, 1) * 100;
end

% arrange pre-determined multiunits
mu = mu(:);
mu = sort(unique([mu, find(spikes.isiViolation > 1)]));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate cluster separation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if all(isfield(spikes, {'mDist', 'iDist', 'lRat'})) && ~force
    warning('separation matrices already calculated. Skipping cluDist');
else
    sprintf('\ncalculating separation matrices\n\n');
        
    % go over each spike group, clean and calc separation
    cnt = 1;
    for igrp = 1 : ngrp
        
        grpid = grps(igrp);

        % load data
        fet = loadNS('datatype', 'fet', 'session', session, 'grpid', grpid);
        clu = loadNS('datatype', 'clu', 'session', session, 'grpid', grpid);

        % find noise and artifact spikes
        rmidx = find(clu <= 1);
        nonClusteredSpks(cnt) = length(rmidx);
        totSpks(cnt) = length(clu);

        % find predetermined mu spikes (mu numbering given by spikes.UID)
        rmclu = spikes.cluID(mu(spikes.shankID(mu) == grpid));
        if ~isempty(rmclu)
            rmidx = unique([rmidx; find(clu == rmclu)]);
        end
        
        uclu = unique(clu);
        nfet = npca * length(spkgrp{grpid});        
        fetDist = [fet, clu];
        fetDist(rmidx, :) = [];
        k = 1;
        for iclu = 1 : length(uclu)
            cluid = uclu(iclu);
            if cluid == 0 || cluid == 1
                continue
            end
            [lRat{cnt}(k, 1), iDist{cnt}(k, 1), ~] =...
                cluDist(fetDist(:, 1 : nfet), find((fetDist(:, end) == cluid)));
            k = k + 1;
        end
        cnt = cnt + 1;
    end
    
    % concatenate cell
    spikes.lRat = cat(1, lRat{:})';
    spikes.iDist = cat(1, iDist{:})';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SU or MU
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
spikes.su = spikes.iDist > 10 & spikes.isiViolation < 1;
spikes.nonClusteredSpks = nonClusteredSpks ./ totSpks;
spikes.suSpks = 1 - length(vertcat(spikes.ts{spikes.su})) / sum(totSpks);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save vars
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if saveVar
    [~, basename] = fileparts(basepath);
    % spikes
    save([basename '.spikes.cellinfo.mat'], 'spikes', '-v7.3')
    
    % cell metrics
    cmName = [basename, '.cell_metrics.cellinfo.mat'];
    if exist(cmName, 'file')
        load(cmName)
        cell_metrics.su = split(num2str(spikes.su))';
        cell_metrics.iDist = spikes.iDist;
        cell_metrics.lRat = spikes.lRat;
        cell_metrics.isiViolation = spikes.isiViolation;       
        save(cmName, 'cell_metrics')
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if graphics
    plotCluster('basepath', basepath, 'spikes', spikes,...
        'vis', vis, 'saveFig', saveFig)
end

end

% EOF