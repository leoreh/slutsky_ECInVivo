function mea = mea_orgNex(varargin)

% gets a .mat file exported from plexon and organized in a struct according
% to cell. assumes fixed order of columns. 
%
% INPUT:
%   fname       char. file name of .mat file. if empty will be extracted
%               from basepath.
%   basepath    char. path to where fname exists {pwd}.
%   forceL      logical. force reload even if mea struct exists {false}
%   saveVar     logical. save mea struct {true}.
%
% DEPENDENCIES
%
% TO DO lIST
%
% 23 dec 21 LH  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addParameter(p, 'fname', '', @ischar);
addParameter(p, 'basepath', pwd, @ischar);
addParameter(p, 'forceL', false, @islogical);
addParameter(p, 'saveVar', true, @islogical);

parse(p, varargin{:})
fname       = p.Results.fname;
basepath    = p.Results.basepath;
forceL      = p.Results.forceL;
saveVar     = p.Results.saveVar;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% indices to nex matrix
idx_ch = 1;
idx_unitID = 2;
idx_spktimes = 3;
wvlength = 29;      % for Anto data wv is 32 samples, and for max 31
fs = 10000;

% file
[~, basename] = fileparts(basepath);
if isempty(fname)
    fname = [basename, '.mat'];
end
meaname = fullfile(basepath, [basename, '.mea.mat']);
if exist(meaname, 'file') && ~forceL
    load(meaname)
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load and organize
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load data
nex = load(fullfile(basepath, fname));

% go through each unit and extract data
unitFields = fieldnames(nex);
nch = length(unitFields);
ncolumns = size(nex.(unitFields{1}), 2);

cnt = 1;
for ich = 1 : nch    
    nunitsCh = unique(nex.(unitFields{ich})(:, idx_unitID));
    
    for iunit = 1 : length(nunitsCh)        
        unitIdx = nex.(unitFields{ich})(:, idx_unitID) == nunitsCh(iunit);      
        mea.ch(cnt) = unique(nex.(unitFields{ich})(unitIdx, idx_ch));
        mea.unitID(cnt) = unique(nex.(unitFields{ich})(unitIdx, idx_unitID));
        mea.spktimes{cnt} = nex.(unitFields{ich})(unitIdx, idx_spktimes);
        wv = nex.(unitFields{ich})(unitIdx, ncolumns - wvlength : ncolumns);
        mea.wv(cnt, :) = mean(wv, 1);
        mea.wv_std(cnt, :) = std(wv, [], 1);
        
        cnt = cnt + 1;
    end
      
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% time bins
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get time bins. mea recordings include only 20 min of every hour. 
binsize = 1200;                                             % 20 min in [s]
lastspike = max(cellfun(@max, mea.spktimes, 'uni', true));     % last spike [s]
tbins = 0 : binsize : lastspike; 
tbins = [tbins, lastspike];

% calculate firing rate in time bins
funh = @(x)  histcounts(x, tbins);
frTemp = cellfun(funh, mea.spktimes, 'uni', false);
fr = [cell2mat(frTemp')] ./ [diff(tbins)];

% manually adjust the timebins. this is because the experimenter sometimes
% changes the number of bins saved per hour. note that tbins_real will
% usually be shorter by 1 element than tbins. to fix this must correct
% spktimes to real times (maybe later). for now this is ok because fr is
% also one element shorter.
bl_tpnt = 4;        % time point when baseline ended [h]
binjumps = 2;       % number of bins skipped once experiment starts
tbins_real = [[1 : bl_tpnt],...
    [5 : binjumps : (length(tbins) - bl_tpnt) * binjumps]];
if length(tbins_real) < length(tbins)
    tbins_real = [tbins_real, tbins_real(end) + binjumps];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mea.fr = fr;
mea.tbins = tbins;
mea.tbins_real = tbins_real * 60 * 60;
mea.info.runtime = datetime(now, 'ConvertFrom', 'datenum');
mea.info.fname = fname;
mea.info.fs = fs;
mea.info.bl_tpnt = bl_tpnt;
mea.info.binjumps = binjumps;

if saveVar
    save(meaname, 'mea')
end

end

% EOF


