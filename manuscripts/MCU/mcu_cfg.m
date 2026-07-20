function cfg = mcu_cfg()

% Cohorts. One list per genotype, aligned with cfg.lbl.grp. mcu_geno maps a
% subject id to its genotype through these lists; a subject in none of them
% is flagged rather than silently treated as a control.
cfg.miceWT = {'lh96', 'lh100', 'lh107', 'lh122', 'lh142', 'lh119', ...
    'lh123', 'lh126'};
cfg.miceMCU = {'lh132', 'lh133', 'lh134', 'lh136', 'lh140', 'lh137'};
cfg.miceCAG = {'raMCU1', 'raMCU2', 'raMCU3', 'raMCU4', 'raMCU5'};

% Colors
cfg.clr.grp = [0.2 0.2 0.2;...      % Control
    0.75 0.55 0.35;...              % MCU-KO (germline)
    0.35 0.55 0.60];                % CAG-MCU-KO (viral)
cfg.clr.bac = [0.5, 0, 0.5, 0.5];
cfg.clr.unit = [0.1, 0.1, 0.4;...   % RS
    0.6, 0.2, 0.2;...               % FS
    0.5, 0.5, 0.5];                 % Other
cfg.clr.cmp = [0.7, 0.2, 0.2;...    % Cyto (red)
    0.2, 0.6, 0.2];                 % Mito (green)

% Variables to load
cfg.vars = {'fr', 'units', 'st_metrics'};

% Variable mapping
cfg.varMap.fr = 'fr.mfr';
cfg.varMap.bRoy = 'st.royer';
cfg.varMap.unitType = 'units.type';

% Labels. Genotype labels must not contain ':' or '_': lme_postHoc splits
% coefficient names on ':' (interaction) and '_' (factor_level), so a label
% carrying either is parsed into a bogus term description.
cfg.lbl.grp = {'Control'; 'MCU-KO'; 'CAG-MCU-KO'};
cfg.lbl.unit = {'RS', 'FS', 'Other'};
cfg.lbl.day = {'BSL'; 'BAC_ON'; 'BAC1'; 'BAC2'; 'BAC3'; 'BAC_OFF'; 'WASH'};

% Save figures path
cfg.savepath = 'D:\OneDrive - Tel-Aviv University\PhD\Slutsky\Manuscripts\Heim2025\Graphs';

end
