
%% ========================================================================
%  LOAD
%  ========================================================================

basepaths = [mcu_basepaths('mea_bac'), mcu_basepaths('mea_mcuko')];
basepaths = [mcu_basepaths('mea_bac')];
nFiles = length(basepaths);

vars = {'mea', 'fr', 'brstDyn', 'brst', 'frRcv', 'frRcv_mdl', ...
    'stats', 'ca', 'prc'};
v = basepaths2vars('basepaths', basepaths, 'vars', vars);

cfg = mcu_cfg;


% -------------------------------------------------------------------------
% % DOWNSAMPLE
% 
% for iFile = 1 : nFiles
% 
%     sig = v(iFile).ca.mito;
%     [nUnits, nSig] = size(sig);
%     nBins = floor(nSig / 60);
%     sig = sig(:, 1 : nBins * 60);
%     sig = reshape(sig, nUnits, 60, nBins);
%     sig = reshape(sum(sig, 2), nUnits, nBins);
%     % sig = fr_denoise(sig, 'frameLen', 30, 'flgPlot', false);
%     v(iFile).ca.mito = sig;
% 
%     sig = v(iFile).ca.cyto;
%     sig = sig(:, 1 : nBins * 60);
%     sig = reshape(sig, nUnits, 60, nBins);
%     sig = reshape(sum(sig, 2), nUnits, nBins);
%     % sig = fr_denoise(sig, 'frameLen', 30, 'flgPlot', false);
%     v(iFile).ca.cyto = sig;
% 
% end

% -------------------------------------------------------------------------
% TABLE

% Vars to align
varMap = struct();
varMap.frt = 'fr.fr';
varMap.btRate = 'brstDyn.rate';
varMap.btDur = 'brstDyn.dur';
varMap.btFreq = 'brstDyn.freq';
varMap.btIBI = 'brstDyn.ibi';
varMap.btFrac = 'brstDyn.bfrac';
% varMap.caCyto = 'ca.cyto';
% varMap.caMito = 'ca.mito';

% Align
% [v, t] = mea_tAlign(v, varMap, 'fr.info.idxPert');
% xVec = t / 3600;

% Unit vars
varMap.uGood = 'fr.uGood';
varMap.uRcv = 'rcv.uRcv';
varMap.uPert = 'rcv.uPert';

% Tag structures
tagFiles.Name = get_mname(basepaths, 0);
tagFiles.Group = repmat(cfg.lbl.grp(1), 1, nFiles);
tagFiles.Group(contains(tagFiles.Name, 'KO')) = cfg.lbl.grp(2);

% Table
tblt = v2tbl('v', v, 'varMap', varMap, 'tagFiles', tagFiles, ...
    'uOffset', 0);

% -------------------------------------------------------------------------
% COMBINE & CLEAN

% Clean bad units
tblt(~tblt.uPert, :) = [];
tblt(~tblt.uGood, :) = [];
tblt.UnitID = categorical(tblt.UnitID);


%% ========================================================================
%  LOG RATIO ANALYSIS
%  ========================================================================

tblu = tblt(:, [1 : 3, 11]);

% Define Windows
idxZero = find(xVec >= 0, 1);
winBsl = [1, idxZero - 8];
winSs = [length(xVec) - winBsl(2), length(xVec) - 1];

% Iterate through Table Variables
tblVars = tblt.Properties.VariableNames;
nTime = length(xVec);

for iVar = 1 : length(tblVars)

    varName = tblVars{iVar};
    val = tblt.(varName);

    % Process only numeric variables that match the time vector length
    if isnumeric(val) && size(val, 2) == nTime

        % Calculate Means within Windows
        valWin = mean(val(:, winBsl(1) : winBsl(2)), 2, 'omitnan');
        valWin(:, 2)  = mean(val(:, winSs(1) : winSs(2)), 2, 'omitnan');
        
        % Add floor to prevent log(0)
        if any(valWin(:) == 0) && all(valWin(:) >= 0)
            c = min(valWin(valWin > 0)) / 2;
            valWin = valWin + c;
            fprintf('[%s] Adding Offset %.4f\n', varName, c);
        end

        % Calculate Log Ratio (Log Fold Change)
        valRcv = (log(valWin(:, 2) ./ valWin(:, 1)));

        % Save as new column in the table
        tblu.([varName, '_bsl']) = valWin(:, 1);
        tblu.([varName, '_ss']) = valWin(:, 2);
        tblu.([varName, '_rcv']) = valRcv;

    end
end

% Remove outliers from mito
% [~, ~, otlBsl] = rmoutliers(tbl.caMito_bsl, "percentiles", [0 98]);
% [~, ~, otlSs] = rmoutliers(tbl.caMito_ss, "percentiles", [0 98]);
% tbl(otlBsl | otlSs, :) = [];

%% ========================================================================
%  PLOT
%  ========================================================================

% Plot Results
% yVar = 'caMito_LR';
% tblGUI_scatHist(tbl, 'grpVar', 'Group', 'yVar', yVar);

% experiments = categories(tbl.Name);
% idxExp = tbl.Name == experiments(2);
% tblGUI_xy(xVec, tbl, 'tileVar', 'Group', 'yVar', 'caMito');



