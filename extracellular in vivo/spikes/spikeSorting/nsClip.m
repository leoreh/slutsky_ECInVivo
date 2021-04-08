function nsClip(varargin)

% creates ns files (spk, res, fet, and clu) from ns files by clipping them.
% clip time is determined according to duration and start time of recording
% (see dur and t below). requires res and spk file in basepath.

% INPUT:
%   basepath    string. path to recording folder {pwd}.
%   spkgrp      array where each cell is the electrodes for a spike group.
%               if empty will be loaded from session info file (cell
%               explorer format)
%   grps        numeric. groups (tetrodes) to work on
%   dur         numeric. duration of trim period [min].
%   t           string. start time of trim period. if empty than will take
%               dur minutes from end of recording (if dur is negative) or
%               dur minutes from start of recording (if dur is positive).
%               can be in the format 'HHmmss' or 'HHmm'.
%   bkup        logical. create backup of files before save in basepath/ns/bkup {false}.
%   fs          numeric. sampling frequency [hz]{20000}
%   nsamps      numeric. number of samples in recording
%
% DEPENDENCIES
%   class2bytes
%
% TO DO LIST
%   # replace memmap w/ fread
%   # clean temp_wh after creation of ns files
%
% 09 nov 20 LH  updates
% 09 mar 20 LH  dur can be pos / neg

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;

p = inputParser;
addOptional(p, 'basepath', pwd);
addOptional(p, 'spkgrp', {}, @iscell);
addOptional(p, 'grps', [], @isnumeric);
addOptional(p, 'dur', [], @isnumeric);
addOptional(p, 't', []);
addOptional(p, 'bkup', false, @islogical);
addOptional(p, 'fs', [], @isnumeric);
addOptional(p, 'nsamps', [], @isnumeric);

parse(p, varargin{:})
basepath    = p.Results.basepath;
spkgrp      = p.Results.spkgrp;
grps        = p.Results.grps;
dur         = p.Results.dur;
t           = p.Results.t;
bkup        = p.Results.bkup;
fs          = p.Results.fs;
nsamps      = p.Results.nsamps;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~, basename] = fileparts(basepath);
sessionfile = fullfile(basepath, [basename '.session.mat']);
if exist(sessionfile, 'file')
    load(sessionfile)
end
if isempty(spkgrp)
    spkgrp = session.extracellular.spikeGroups.channels;
end
if isempty(fs)
    fs = session.extracellular.sr;
end
if isempty(nsamps)
    nsamps = session.general.nsamps;
end
if isempty(grps)
    grps = 1 : length(spkgrp);
end
ngrps = length(grps);

% clip clu only if exists
if isempty(dir('*clu*'))
    cluFlag = false;
else
    cluFlag = true;
end

recStart = '';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find trim boundries [samples]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dur = dur * 60 * fs;
if ischar(t)
    dt = tstamp2time('dtstr', t, 'fs', fs);
    if ~isempty(dur)
        trimEdges = sort([dt dt + dur]);
    else
        trimEdges = [dt nsamps];
    end
else
    if dur < 0
        trimEdges = [nsamps + dur nsamps];
    elseif dur > 0
        trimEdges = [1 dur];
    else
        trimEdges = [1 nsamps];
    end
end
if trimEdges(1) < 1
    warning('Requested time beyond recording duration')
    trimEdges(1) = 1;
elseif trimEdges(2) > nsamps
    warning('Requested time beyond recording duration')
    trimEdges(2) = nsamps;
end

% save trim info
infoname = fullfile(basepath, [basename, '.datInfo.mat']);
if exist(infoname, 'file')
    load(infoname)
end
datInfo.spktrim.dur = dur / fs / 60;
datInfo.spktrim.edges = trimEdges;
datInfo.spktrim.recStart = recStart;
datInfo.spktrim.t = t;
save(infoname, 'datInfo', '-v7.3');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% go over groups, trim and save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for igrp = 1 : ngrps
    grp = grps(igrp);
    grpchans = spkgrp{igrp};
    
    % load data
    if cluFlag
        clu = loadNS('datatype', 'clu', 'session', session, 'grpid', grp);
    end
    res = loadNS('datatype', 'res', 'session', session, 'grpid', grp);
    nspks = length(res);
    spk = loadNS('datatype', 'spk', 'session', session, 'grpid', grp,...
        'nspks', nspks);
    
    % find idx of spikes in trim boundries
    spkidx = find(res > trimEdges(1) & res < trimEdges(2));
    nspks = length(spkidx);
    
    % clip data
    res = res(spkidx);
    if cluFlag; clu = clu(spkidx); end
    spk = spk(:, :, spkidx);
    
    % recalculate fet
    fprintf('Computing PCAs...')
    nfet = length(grpchans) * 3 + length(grpchans) + 1;
    fetMat = zeros(nspks, nfet);
    enrgIdx = length(grpchans) * 3;
    if ~isempty(spk)
        for ichan = 1 : length(grpchans)
            [~, pcFeat] = pca(permute(spk(ichan, :, :), [3, 2, 1]), 'NumComponents', 3);
            chEnrgy = sum(abs(permute(spk(ichan, :, :), [3, 2, 1])), 2);
            fetMat(:, ichan * 3 - 2 : ichan * 3) = pcFeat;
            fetMat(:, enrgIdx + ichan) = (chEnrgy);
        end
    end
    fetMat(:, end) = double(res);
    fet = int32(fetMat');
    fprintf('done. \n')
    
    % save data
    if cluFlag
        saveNS(clu, 'datatype', 'clu', 'session', session, 'grpid', grp, 'bkup', bkup);
    end
    saveNS(res, 'datatype', 'res', 'session', session, 'grpid', grp, 'bkup', bkup);
    saveNS(spk, 'datatype', 'spk', 'session', session, 'grpid', grp, 'bkup', bkup);
    saveNS(fet, 'datatype', 'fet', 'session', session, 'grpid', grp, 'bkup', bkup);
    
end

end

% EOF