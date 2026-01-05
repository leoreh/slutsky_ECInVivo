function mea = mea_orgNex(varargin)

% gets a .mat file exported from plexon and organized in a struct according
% to cell. assumes fixed order of columns. 
%
% INPUT:
%   fname       char. file name of .mat file. if empty will be extracted
%               from basepath.
%   basepath    char. path to where fname exists {pwd}.
%   flgForce    logical. force reload even if mea struct exists {false}
%   flgSave     logical. save mea struct {true}.
%   flgPlot     logical. plot fr and waveforms {false}.
%
% 23 dec 21 LH  

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addParameter(p, 'fname', '', @ischar);
addParameter(p, 'basepath', pwd, @ischar);
addParameter(p, 'flgForce', false, @islogical);
addParameter(p, 'flgSave', true, @islogical);
addParameter(p, 'flgPlot', false, @islogical);

parse(p, varargin{:})
fname       = p.Results.fname;
basepath    = p.Results.basepath;
flgForce    = p.Results.flgForce;
flgSave     = p.Results.flgSave;
flgPlot     = p.Results.flgPlot;

%% ========================================================================
%  PREPERATIONS
%  ========================================================================

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
if exist(meaname, 'file') && ~flgForce
    load(meaname)
    return
end

%% ========================================================================
%  LOAD
%  ========================================================================

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

%% ========================================================================
%  SAVE
%  ========================================================================

mea.info.runtime = datetime("now");
mea.info.fname = fname;
mea.info.fs = fs;

if flgSave
    save(meaname, 'mea')
end

%% ========================================================================
%  PLOT
%  ========================================================================

if flgPlot

    % Time bins
    binsize = 1200;                                                % 20 min in [s]
    lastspike = max(cellfun(@max, mea.spktimes, 'uni', true));
    tbins = 0 : binsize : lastspike;
    tbins(end) = lastspike;
    tVec = tbins(1:end-1) + diff(tbins)/2;

    % calculate firing rate in time bins
    funh = @(x)  histcounts(x, tbins);
    frTemp = cellfun(funh, mea.spktimes, 'uni', false);
    fr = [cell2padmat(frTemp, 1)] ./ [diff(tbins)];

    % create table
    tbl = table();
    tbl.unitID = categorical(1 : size(fr, 1))';
    tbl.spktimes = mea.spktimes';
    tbl.fr = fr;
    tbl.wv = mea.wv;

    % plot fr
    tblGUI_xy(tVec, tbl, 'yVar', 'fr', 'grpVar', 'unitID', 'xLbl', 'Time (s)');

    % plot wv
    tWv = (0:size(mea.wv, 2)-1) / mea.info.fs * 1000;
    tblGUI_xy(tWv, tbl, 'yVar', 'wv', 'xLbl', 'Time (ms)');

    % plot raster
    tblGUI_raster(tbl, 'timesVar', 'spktimes');
    
end

end     % EOF
