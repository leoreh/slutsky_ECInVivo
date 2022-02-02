% mea_sessions

% loads all relavent files from multiple sessions (experiments) and does
% stuff

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get all folders in masterpath
masterpath = 'G:\Data\MEA\max';
d = dir(masterpath);
d = d([d(:).isdir]);
d = d(~ismember({d(:).name},{'.','..'}));
basenames = {d.name};
for isession = 1 : length(basenames)
    basepaths{isession} = fullfile(masterpath, basenames{isession});
end

% get specific basepaths
basepaths = {'D:\Google Drive\Data\mea\baclofen\210325_082300',...
    'D:\Google Drive\Data\mea\baclofen\210527_045600',...
    'D:\Google Drive\Data\mea\baclofen\190824_123000',...
    'D:\Google Drive\Data\mea\baclofen\190910_045200',...
    'D:\Google Drive\Data\mea\ketamine\210516_174646',...
    'D:\Google Drive\Data\mea\ketamine\201227_161755',...
        };
       
nsessions = length(basepaths);

% analyze all sessions
for isession = 1 : nsessions
    %     mea_analyze('basepath', basepaths{isession},...
    %         'winBL', [0, 120 * 60], 'graphics', true, 'forceA', false)
    
    winCalc = v(isession).monosyn.info.winCalc;
    mea = v(isession).mea;
    monosyn = monoSyn_wrapper('spktimes', mea.spktimes, 'basepath', pwd,...
        'winCalc', winCalc, 'saveVar', true, 'graphics', false,...
        'forceA', true, 'fs', mea.info.fs, 'saveFig', false,...
        'wv', mea.wv, 'wv_std', mea.wv_std);    
    
end

% load vars from each session
varsFile = ["fr"; "mea"; "st_metrics"; "swv_metrics"; "cell_metrics"; "monoSyn.mat"];
varsName = ["fr"; "mea"; "st"; "swv"; "cm"; "monosyn"];
v = getSessionVars('basepaths', basepaths, 'varsFile', varsFile,...
    'varsName', varsName);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% concat data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize
tp = []; spkw = []; royer = []; lidor = []; mfr =[]; tau_rise = [];
mizuseki = []; lvr = []; asym = []; hpk = []; rs = []; fs = [];
slopeTail = []; eiDom = []; slopeTp = []; ampTail = [];


for isession = 1 : nsessions
    
    mfr = [mfr, v(isession).fr.mfr'];
    
    eiDom = [eiDom; v(isession).monosyn.eiDom];
    
    asym = [asym, v(isession).swv.asym];
    hpk = [hpk, v(isession).swv.hpk];
    tp = [tp, v(isession).swv.tp];
    spkw = [spkw, v(isession).swv.spkw];
    slopeTail = [slopeTail, v(isession).swv.slopeTail];
    slopeTp = [slopeTp, v(isession).swv.slopeTp];
    ampTail = [ampTail, v(isession).swv.ampTail];
    
    lvr = [lvr, v(isession).st.lvr];
    royer = [royer, v(isession).st.royer];
    lidor = [lidor, v(isession).st.lidor];
    mizuseki = [mizuseki, v(isession).st.mizuseki];
    tau_rise = [tau_rise, v(isession).st.tau_rise];
    
end

% mfr = normalize(mfr, 'range', [0.1 1]);

% separate units according to E-I dominance
eiThr = [-0.3, 0.3];
eu = eiDom > eiThr(2);
iu = eiDom < eiThr(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ---------------------------------------------------------------------
% classification

fh = figure;
subplot(2, 3, 1)
plot(mfr, eiDom, '.k', 'markerSize', 15)
hold on
plot(mfr(eu), eiDom(eu), '.b', 'markerSize', 15)
plot(mfr(iu), eiDom(iu), '.r', 'markerSize', 15)
set(gca, 'XScale', 'log')
plot(xlim, [eiThr(1) eiThr(1)], '--r')
plot(xlim, [eiThr(2) eiThr(2)], '--b')
xlabel('Firing Rate')
ylabel('E-I dominance')

subplot(2, 3, 2)
var1 = tp;
var2 = royer;
plot(var1, var2, '.k', 'markerSize', 15)
hold on
plot(var1(eu), var2(eu), '.b', 'markerSize', 15)
plot(var1(iu), var2(iu), '.r', 'markerSize', 15)
set(gca, 'YScale', 'log')
xlabel('TP')
ylabel('Royer')

subplot(2, 3, 3)
var1 = lvr;
var2 = spkw;
plot(var1, var2, '.k', 'markerSize', 15)
hold on
plot(var1(eu), var2(eu), '.b', 'markerSize', 15)
plot(var1(iu), var2(iu), '.r', 'markerSize', 15)
xlabel('LvR')
ylabel('Spike Width')

subplot(2, 3, 4)
var1 = mfr;
var2 = lvr;
plot(var1, var2, '.k', 'markerSize', 15)
hold on
plot(var1(eu), var2(eu), '.b', 'markerSize', 15)
plot(var1(iu), var2(iu), '.r', 'markerSize', 15)
xlabel('MFR')
ylabel('LvR')
set(gca, 'XScale', 'log')

subplot(2, 3, 5)
var1 = tp;
var2 = asym;
plot(var1, var2, '.k', 'markerSize', 15)
hold on
plot(var1(eu), var2(eu), '.b', 'markerSize', 15)
plot(var1(iu), var2(iu), '.r', 'markerSize', 15)
xlabel('TP')
ylabel('Asymmetry')

subplot(2, 3, 6)
var1 = slopeTp;
var2 = ampTail;
plot(var1, var2, '.k', 'markerSize', 15)
hold on
plot(var1(eu), var2(eu), '.b', 'markerSize', 15)
plot(var1(iu), var2(iu), '.r', 'markerSize', 15)
xlabel('slopeTp')
ylabel('ampTail')

% 
% % save
% figname = fullfile(masterpath, 'cellClass');
% export_fig(figname, '-jpg', '-transparent', '-r300')

