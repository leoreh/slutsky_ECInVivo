%% ========================================================================
%  MCU BEHAVIOR ANALYSIS
%  ========================================================================
% Recognition memory in CA1:MCU-KO and CA2:MCU-KO mice across two tasks:
% Novel Object Recognition (NOR) and Social Recognition. Computes
% per-animal discrimination indices, effect sizes, the subfield x genotype
% interaction (the formal dissociation test), and equivalence testing on
% the null arm (CA2 NOR).
%
% Source data: BehaviorData.xlsx, sheet 'Long' (long-format table built
% from the four per-task raw sheets).
%
% Output: mcu_behaviorStats.xlsx in the manuscript Results folder, with
% one sheet per result block (DI, EffSize, Interaction, Equivalence,
% Power).
%
% Dependencies: effSize_d, effSize_tost, effSize_pwr (utilities/),
% guiTbl_bar (graphics/).


%% ========================================================================
%  CONFIG
%  ========================================================================

% Paths
behPath  = 'D:\OneDrive - Tel-Aviv University\PhD\Slutsky\Manuscripts\MCU\Results\Data_ShayA';
behXls   = 'BehaviorData.xlsx';
outPath  = 'D:\OneDrive - Tel-Aviv University\PhD\Slutsky\Manuscripts\MCU\Results';
outXls   = 'mcu_behaviorStats.xlsx';

% Flags
flgPlot  = false;
flgSave  = true;

% Statistical settings
alpha    = 0.05;
confLvl  = 0.95;
sesoiD   = 0.5;             % equivalence bound, Cohen's d units


%% ========================================================================
%  LOAD
%  ========================================================================

tblBeh = readtable(fullfile(behPath, behXls), 'Sheet', 'Long', ...
    'TextType', 'string');

% Categoricals (ordered for sensible plotting)
tblBeh.Subfield = categorical(tblBeh.Subfield, {'CA1', 'CA2'});
tblBeh.Genotype = categorical(tblBeh.Genotype, {'Control', 'MCU-KO'});
tblBeh.Task     = categorical(tblBeh.Task,     {'Object', 'Social'});
tblBeh.Phase    = categorical(tblBeh.Phase);

% Exclude column
exAnimals = unique(tblBeh.sbjID(tblBeh.Exclude == 1));
tblBeh    = tblBeh(~ismember(tblBeh.sbjID, exAnimals), :);
fprintf('Excluded %d animals: %s\n', numel(exAnimals), strjoin(exAnimals, ', '));


%% ========================================================================
%  DERIVED MEASURES (per animal)
%  ========================================================================
% Object Recognition (Test phase):
%   DI    = (T_test_B - T_test_A) / (T_test_B + T_test_A)
% Social Recognition:
%   DisHI = (T_B1 - T_A2) / (T_B1 + T_A2)     Dishabituation Index
%   HI    = (T_A1 - T_A2) / (T_A1 + T_A2)     Habituation Index
% DI is populated for Object only; DisHI/HI for Social only. The
% dissociation test compares DI (Object) and HI (Social) across
% subfields and genotypes, matching the primary index reported in
% Figure 3F and the Methods. DisHI is preserved alongside HI for
% transparency and any auxiliary analysis.

sbjIDs = unique(tblBeh.sbjID);
nSbj   = numel(sbjIDs);
varNms = {'sbjID', 'Subfield', 'Genotype', 'Task', ...
          'DI', 'DisHI', 'HI', 'T_familiar', 'T_novel', 'T_total'};
varTyp = {'string', 'categorical', 'categorical', 'categorical', ...
          'double', 'double', 'double', 'double', 'double', 'double'};
tblDI  = table('Size', [nSbj, numel(varNms)], 'VariableTypes', varTyp, ...
    'VariableNames', varNms);

for iS = 1:nSbj
    rows  = tblBeh(tblBeh.sbjID == sbjIDs(iS), :);
    sub   = rows.Subfield(1);
    geno  = rows.Genotype(1);
    task  = rows.Task(1);

    if task == 'Object'
        tFam = rows.Time(rows.Phase == 'Test' & rows.Target == "A");
        tNov = rows.Time(rows.Phase == 'Test' & rows.Target == "B");
        tTot = tFam + tNov;
        di    = (tNov - tFam) / tTot;
        disHI = NaN;
        hi    = NaN;

    else                                            % Social
        tA1 = rows.Time(rows.Target == "A1");
        tA2 = rows.Time(rows.Target == "A2");
        tB1 = rows.Time(rows.Target == "B1");
        tFam = tA2;
        tNov = tB1;
        tTot = tA1 + tA2 + tB1;
        di    = NaN;
        disHI = (tB1 - tA2) / (tB1 + tA2);
        hi    = (tA1 - tA2) / (tA1 + tA2);
    end

    tblDI(iS, :) = {sbjIDs(iS), sub, geno, task, di, disHI, hi, tFam, tNov, tTot};
end

% Drop rows where the task-relevant primary index is missing
keep = (tblDI.Task == 'Object' & ~isnan(tblDI.DI)) | ...
       (tblDI.Task == 'Social' & ~isnan(tblDI.DisHI));
tblDI = tblDI(keep, :);


%% ========================================================================
%  SANITY PLOTS (guiTbl_bar)
%  ========================================================================

if flgPlot
    % Object DI x Subfield, grouped by Genotype
    guiTbl_bar(tblDI(tblDI.Task == 'Object', :), ...
        'yVar', 'DI', 'xVar', 'Subfield', 'grpVar', 'Genotype');
    set(gcf, 'Name', 'DI - Object');

    % Social Habituation Index x Subfield, grouped by Genotype
    guiTbl_bar(tblDI(tblDI.Task == 'Social', :), ...
        'yVar', 'HI', 'xVar', 'Subfield', 'grpVar', 'Genotype');
    set(gcf, 'Name', 'HI - Social');

    % Total exploration time (gross-engagement check; flags locomotor or
    % motivational confounds before they enter the index interpretation)
    guiTbl_bar(tblDI, 'yVar', 'T_total', 'xVar', 'Subfield', 'grpVar', 'Genotype');
    set(gcf, 'Name', 'Total exploration');
end


%% ========================================================================
%  EFFECT SIZES (Cohen's d, Hedges' g, 95% CI)
%  ========================================================================
% Compute Ctrl-vs-KO effect size for each (subfield, task) cell.

combos = {'CA1', 'Object'; 'CA1', 'Social'; ...
          'CA2', 'Object'; 'CA2', 'Social'};
nC = size(combos, 1);

effNms = {'Subfield', 'Task', 'n_Ctrl', 'n_KO', ...
          'mean_Ctrl', 'mean_KO', 'Cohen_d', 'Hedges_g', ...
          'd_CI_low', 'd_CI_high'};
effTyp = {'categorical', 'categorical', 'double', 'double', ...
          'double', 'double', 'double', 'double', 'double', 'double'};
tblEff = table('Size', [nC, numel(effNms)], 'VariableTypes', effTyp, ...
    'VariableNames', effNms);

for iC = 1:nC
    sub  = combos{iC, 1};
    task = combos{iC, 2};
    if strcmp(task, 'Object'), yCol = 'DI'; else, yCol = 'HI'; end
    vC = tblDI.(yCol)(tblDI.Subfield == sub & tblDI.Task == task & tblDI.Genotype == 'Control');
    vK = tblDI.(yCol)(tblDI.Subfield == sub & tblDI.Task == task & tblDI.Genotype == 'MCU-KO');

    [d, g, ci, ~] = effSize_d(vC, vK, 'confLvl', confLvl);

    tblEff(iC, :) = {sub, task, numel(vC), numel(vK), ...
        mean(vC), mean(vK), d, g, ci(1), ci(2)};
end

disp('Effect sizes (Ctrl vs KO):'); disp(tblEff);


%% ========================================================================
%  INTERACTION (subfield x genotype, per task)
%  ========================================================================
% Each animal contributes one DI per task; this is a 2x2 between-subjects
% design. Use OLS via fitlm. The subfield:genotype interaction is the
% formal statement of the dissociation.

% Object
tblObj = tblDI(tblDI.Task == 'Object', :);
mdlObj = fitlm(tblObj, 'DI ~ Subfield * Genotype');
anvObj = anova(mdlObj);
anvObj.Term = string(anvObj.Properties.RowNames);
anvObj.Properties.RowNames = {};
anvObj = movevars(anvObj, 'Term', 'Before', 1);

% Social
tblSoc = tblDI(tblDI.Task == 'Social', :);
mdlSoc = fitlm(tblSoc, 'HI ~ Subfield * Genotype');
anvSoc = anova(mdlSoc);
anvSoc.Term = string(anvSoc.Properties.RowNames);
anvSoc.Properties.RowNames = {};
anvSoc = movevars(anvSoc, 'Term', 'Before', 1);

% Combined ANOVA table with task label
anvObj.Task = repmat(categorical("Object"), height(anvObj), 1);
anvSoc.Task = repmat(categorical("Social"), height(anvSoc), 1);
tblAnv = [anvObj; anvSoc];
tblAnv = movevars(tblAnv, 'Task', 'Before', 1);

disp('Interaction ANOVA (per task):'); disp(tblAnv);

% Pull interaction p-values for headline reporting
pIntObj = anvObj.pValue(strcmp(anvObj.Term, 'Subfield:Genotype'));
pIntSoc = anvSoc.pValue(strcmp(anvSoc.Term, 'Subfield:Genotype'));
fprintf('Subfield x Genotype interaction:  Object p = %.4f  |  Social p = %.4f\n', ...
    pIntObj, pIntSoc);


%% ========================================================================
%  EQUIVALENCE (TOST on CA2 NOR)
%  ========================================================================
% The dissociation needs the CA2 NOR null to be defended, not just
% non-significant. SESOI = max(sesoiD, 0.5 * |CA1 NOR d|).

dCA1obj = tblEff.Cohen_d(tblEff.Subfield == 'CA1' & tblEff.Task == 'Object');
sesoiUse = max(sesoiD, 0.5 * abs(dCA1obj));

ctrlCA2obj = tblDI.DI(tblDI.Subfield == 'CA2' & tblDI.Task == 'Object' & tblDI.Genotype == 'Control');
koCA2obj   = tblDI.DI(tblDI.Subfield == 'CA2' & tblDI.Task == 'Object' & tblDI.Genotype == 'MCU-KO');
[pTost, flgEq, infoTost] = effSize_tost(ctrlCA2obj, koCA2obj, sesoiUse, ...
    'unit', 'd', 'alpha', alpha);

tblTost = table(categorical("CA2"), categorical("Object"), ...
    sesoiUse, infoTost.diff, infoTost.se_diff, infoTost.ci90(1), infoTost.ci90(2), ...
    infoTost.p_lower, infoTost.p_upper, pTost, flgEq, ...
    'VariableNames', {'Subfield', 'Task', ...
        'SESOI_d', 'diff_raw', 'se_diff', 'ci90_low', 'ci90_high', ...
        'p_lower', 'p_upper', 'p_TOST', 'Equivalent'});

fprintf('TOST (CA2 NOR, SESOI = %.2f d):  p_TOST = %.4f  -> equivalent = %s\n', ...
    sesoiUse, pTost, mat2str(flgEq));


%% ========================================================================
%  POWER (detectable d at current n, per cohort)
%  ========================================================================
% For each (subfield, task), report the minimum detectable Cohen's d at
% 80% power and current sample size. Useful in legends to qualify nulls.

pwrTarget = 0.80;
nPwr      = height(tblEff);
detD      = nan(nPwr, 1);
for iC = 1:nPwr
    nVec = [tblEff.n_Ctrl(iC), tblEff.n_KO(iC)];
    % Bisection on d for the smallest effect detectable at pwrTarget,
    % using the exact unequal-n formula.
    lo = 0.01; hi = 5;
    for it = 1:60
        mid = 0.5 * (lo + hi);
        [~, info] = effSize_pwr(mid, 'n', nVec, 'alpha', alpha);
        if info.pwr < pwrTarget, lo = mid; else, hi = mid; end
        if hi - lo < 1e-4, break; end
    end
    detD(iC) = hi;
end
tblPwr = tblEff(:, {'Subfield', 'Task', 'n_Ctrl', 'n_KO'});
tblPwr.power_target  = repmat(pwrTarget, nPwr, 1);
tblPwr.alpha         = repmat(alpha, nPwr, 1);
tblPwr.d_detectable  = detD;

disp('Detectable d at current n (80% power):'); disp(tblPwr);


%% ========================================================================
%  EXPORT
%  ========================================================================

if flgSave
    fullOut = fullfile(outPath, outXls);
    if isfile(fullOut), delete(fullOut); end

    writetable(tblDI,   fullOut, 'Sheet', 'DI');
    writetable(tblEff,  fullOut, 'Sheet', 'EffSize');
    writetable(tblAnv,  fullOut, 'Sheet', 'Interaction');
    writetable(tblTost, fullOut, 'Sheet', 'Equivalence');
    writetable(tblPwr,  fullOut, 'Sheet', 'Power');

    fprintf('Wrote results to %s\n', fullOut);
end
