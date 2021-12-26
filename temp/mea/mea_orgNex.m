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
addOptional(p, 'fname', '', @ischar);
addOptional(p, 'basepath', pwd, @ischar);
addOptional(p, 'forceL', false, @islogical);
addOptional(p, 'saveVar', true, @islogical);

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
wvlength = 30;
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
    
    for iunit = 1 : nunitsCh
        
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
% save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mea.info.runtime = datetime(now, 'ConvertFrom', 'datenum');
mea.info.fname = fname;
mea.info.fs = fs;

if saveVar
    save(meaname, 'mea')
end

end

% EOF


