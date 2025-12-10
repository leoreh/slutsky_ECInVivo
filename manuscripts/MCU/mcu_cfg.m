function cfg = mcu_cfg()

cfg.clr.grp = [0.2 0.2 0.2;...
    0.75 0.55 0.35];
cfg.clr.bac = [0.5, 0, 0.5, 0.5];
cfg.clr.unitType = [0.1, 0.1, 0.4;...   % RS
    0.6, 0.2, 0.2;...               % FS
    0.5, 0.5, 0.5];                 % Other

% Variable mapping
cfg.varMap.FR = 'fr.mfr';
cfg.varMap.BRoy = 'st.royer';
cfg.varMap.BLidor = 'st.lidor';
% cfg.varMap.BMiz = 'st.mizuseki';
% cfg.varMap.BSpks = 'brst.bspks';
cfg.varMap.UnitType = 'units.type';

% Variables to load
% cfg.vars = {'fr', 'units', 'st_metrics', 'st_brst'};
cfg.vars = {'fr', 'units', 'st_metrics'};

% Labels
cfg.lbl.grp = {'Control'; 'MCU-KO'};
cfg.lbl.unit = {'RS', 'FS', 'Other'};
cfg.lbl.day = {'BSL'; 'BAC1'; 'BAC2'; 'BAC3'; 'WASH'};

% remove bac on and off
cfg.idxRm = [2, 6];

% Save figures path
cfg.savepath = 'D:\OneDrive - Tel-Aviv University\PhD\Slutsky\Manuscripts\Heim2025\Graphs';

end
