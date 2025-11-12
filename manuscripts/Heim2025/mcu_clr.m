function [clr, lbl] = mcu_clr()

clr.grp = [0.2 0.2 0.2;...
    0.75 0.55 0.35];
clr.bac = [0.5, 0, 0.5, 0.5];
clr.unitType = [0.1, 0.1, 0.4;...   % pPYR
    0.6, 0.2, 0.2];                 % pINT

% Variable mapping
lbl.varMap.FR = 'fr.mfr';
lbl.varMap.BRoy = 'st.royer';
lbl.varMap.BLidor = 'st.lidor';
lbl.varMap.BMiz = 'st.mizuseki';
lbl.varMap.BSpks = 'brst.bspks';
lbl.varMap.UnitType = 'units.clean';
lbl.varMap.PRC = 'prc.prc0_norm';

% Variables to load
lbl.vars = {'fr', 'units', 'st_metrics', 'st_brst', 'prc'};

% Labels
lbl.grp = {'Control'; 'MCU-KO'};
lbl.unit = {'RS', 'FS'};
lbl.day = {'BSL'; 'BAC1'; 'BAC2'; 'BAC3'; 'WASH'};

% remove bac on and off
lbl.idxRm = [2, 6];

% Save figures path
lbl.savepath = 'D:\OneDrive - Tel-Aviv University\PhD\Slutsky\Manuscripts\Heim2025\Graphs';

end
