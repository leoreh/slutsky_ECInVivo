
% LOAD
load('tbl.mat')
load('xVec.mat')


% PLOTS
guiVars = {'Name', 'Group', 'UnitID', 'frt', 'btRate', 'btDur', 'btFreq', 'btIBI', 'btFrac'};
guiTbl_xy(xVec, tbl(:, guiVars));

guiVars = {'Name', 'Group', 'UnitID', 'bRate', 'bDur', 'bFreq', 'bIBI', 'bFrac', 'fr', 'frSs', 'spkDfct', 'rcvTime'};
guiTbl_scatHist(tbl(:, guiVars), 'xVar', 'fr', 'yVar', 'bFrac', 'grpVar', 'Group');

guiTbl_bar(tbl(:, guiVars), 'yVar', 'bFrac', 'xVar', 'Group');

guiVars = {'Name', 'Group', 'UnitID', 'spktimes'};
guiTbl_raster(tbl(:, guiVars), 'grpVar', 'Name', 'grpVal', 'ctrl1')