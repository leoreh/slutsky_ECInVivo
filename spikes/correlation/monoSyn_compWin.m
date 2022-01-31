function monosyn = monoSyn_compWin(varargin)

% compare STGs from two time windows 
%
% INPUT:
%   spktimes    cell of spktimes [s].  
%   basepath    char. path to recording session.
%   winCalc     2-element vector describing the window for calculating the
%               stgs [s] {[0 Inf]}.
%   fs          numeric. sampling frequency {10000}.
%   wv          2 x n numeric. mean waveform of units.
%   wv_std      2 x n numeric. std of units waveforms.
%   saveVar     logical / char. save variable {true}. if char then save
%               name will include saveVar as a suffix
%   graphics    logical. plot graphics or not {true}. 
%   forceA      logical. reanalyze even if struct file exists
%
% DEPENDENCIES
%   monoSyn_wrapper    
%
% TO DO lIST
%
% 31 jan 22 LH  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% validate_win = @(win) assert(isnumeric(win) && length(win) == 2,...
%     'time window must be in the format [start end]');
% 
% p = inputParser;
% addOptional(p, 'spktimes', {}, @iscell);
% addOptional(p, 'basepath', pwd, @ischar);
% addOptional(p, 'winCalc', [0 Inf], validate_win);
% addOptional(p, 'fs', 10000, @isnumeric);
% addOptional(p, 'wv', [], @isnumeric);
% addOptional(p, 'wv_std', [], @isnumeric);
% addOptional(p, 'saveVar', true);
% addOptional(p, 'graphics', true, @islogical);
% addOptional(p, 'forceA', false, @islogical);
% 
% parse(p, varargin{:})
% spktimes    = p.Results.spktimes;
% basepath    = p.Results.basepath;
% winCalc     = p.Results.winCalc;
% fs          = p.Results.fs;
% wv          = p.Results.wv;
% wv_std      = p.Results.wv_std;
% saveVar     = p.Results.saveVar;
% graphics    = p.Results.graphics;
% forceA      = p.Results.forceA;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% session params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

basepath = pwd;
cd(basepath)
[~, basename] = fileparts(basepath);

frfile = fullfile(basepath, [basename, '.fr.mat']);
meafile = fullfile(basepath, [basename, '.mea.mat']);
load(frfile)
load(meafile)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% analyze stg in time windows
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plot of MFR 
fh = figure;
yLimit = [0 ceil(max(mean(fr.strd, 2, 'omitnan')))];
hold on
plot(fr.tstamps / 60 / 60, mean(fr.strd, 1, 'omitnan'), 'k', 'LineWidth', 2)
axis tight
xlabel('Time [~h]')
ylabel('Firing rate [Hz]')
set(gca, 'box', 'off')

% time windows
winD = 60 * 60;                 % winCalc size [s]
winStart = [0.2, 7.5]';     % start time for each win [h]
winCalc = [winStart * 60 * 60, winStart * 60 * 60 + winD];
nwin = size(winCalc, 1);

% mono synaptic interactions (spike transmission gain)
for iwin = 1 : nwin
    winSave = round(winCalc(iwin, :) / 60 / 60);
    saveVar = [num2str(winSave(1)), '-', num2str(winSave(2))];
    ms(iwin) = monoSyn_wrapper('spktimes', mea.spktimes, 'basepath', pwd,...
        'winCalc', winCalc(iwin, :), 'saveVar', saveVar, 'graphics', false,...
        'forceA', true, 'fs', mea.info.fs, 'wv', mea.wv, 'wv_std', mea.wv_std);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load ms structs if already analzed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

monofiles = dir('*monoSyn*');
monofiles = natsort({monofiles.name});
clear ms
for ifile = 1 : length(monofiles)
    load(monofiles{ifile})
    ms(ifile) = monosyn;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compare
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nunits = length(ms(1).eiDom);

% initialize
eStg = cell(1, nwin);
iStg = cell(1, nwin);
nspks = [];
for iwin = 1 : nwin
    
    % significant synapses from win
    [u1, u2] = find(ms(iwin).eSig);
    uE = sub2ind([nunits, nunits], u1, u2);
    [u1, u2] = find(ms(iwin).iSig);
    uI = sub2ind([nunits, nunits], u1, u2);
    
    for iwin2 = 1 : nwin
        
        eStg{iwin} = [eStg{iwin}, ms(iwin2).eStg(uE)];
        iStg{iwin} = [iStg{iwin}, ms(iwin2).iStg(uI)];
        
    end
    
    % nspks per cell
    nspks = [nspks, ms(iwin).nspks'];
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compare
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clr = 'rgb';
clear lgndTxt
fh = figure;

subplot(2, 2, [1, 2])
plot(fr.tstamps / 60 / 60, mean(fr.strd, 1, 'omitnan'), 'k', 'LineWidth', 2)
hold on
for iwin = 1 : nwin
    plot([ms(iwin).info.winCalc; ms(iwin).info.winCalc] / 60 / 60,...
        ylim, '--', 'Color', clr(iwin), 'LineWidth', 2)
end
axis tight
xlabel('Time [~h]')
ylabel('Firing rate [Hz]')
set(gca, 'box', 'off')
title(basename)

subplot(2, 2, 3)
hold on
for iwin = 1 : nwin
    plot([1 : nwin], mean(eStg{iwin}, 1, 'omitnan'),...
        'Color', clr(iwin))
    title('excitatory')
    ylim([0, 0.1])
    lgndTxt{iwin} = sprintf('ref%d = %d synapses', iwin, length(eStg{iwin}));
end
legend(lgndTxt, 'location', 'northwest')
ylabel('E-STG')
xticks([1 : nwin])

subplot(2, 2, 4)
hold on
for iwin = 1 : nwin
    plot([1 : nwin], mean(iStg{iwin}, 1, 'omitnan'),...
        'Color', clr(iwin))
    title('inhibitory')
    ylim([-0.01, 0])
    lgndTxt{iwin} = sprintf('ref%d = %d synapses', iwin, length(iStg{iwin}));
end
legend(lgndTxt, 'location', 'southeast')
ylabel('I-STG')





