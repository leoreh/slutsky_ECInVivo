% fr_sessions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
forceL = false;
forceA = false;

% should allow user to input varName or columnn index
colName = 'Session';                    % column name in xls sheet where dirnames exist
vars = ["cell_metrics.cellinfo";...
    "spikes.cellinfo"];      % string array of variables to load
cond = ["manCur"];                      % column name of logical values for each session. only if true than session will be loaded. can be a string array and than all conditions must be met.
basepath = 'H:\Data\Processed\lh52';
sessionlist = 'sessionList.xlsx';       % must include extension
fs = 20000;                             % can also be loaded from datInfo

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get directory paths
if exist('dirnames', 'var') && isstring(dirnames)
    % ALT 1: user input dirnames
    dirnames = dirnames;
elseif ischar(sessionlist) && contains(sessionlist, 'xlsx')
    % ALT 2: get dirnames from xlsx file
    sessionInfo = readtable(fullfile(basepath, sessionlist));
    icol = strcmp(sessionInfo.Properties.VariableNames, colName);
    dirnames = string(table2cell(sessionInfo(:, icol)));
    % check dirnames meet conditions
    for i = 1 : length(cond)
        iicol(i) = find(strcmp(sessionInfo.Properties.VariableNames, char(cond(i))));
        irow{i} = find(sessionInfo{:, iicol(i)} == 1);
    end
    if length(irow) > 1
        dirnames = dirnames(intersect(irow{:}));
    else
        dirnames = dirnames(irow{1});
    end
end

nsessions = length(dirnames);

% load files
if forceL   
    d = cell(length(dirnames), length(vars));
    for i = 1 : nsessions
        filepath = char(fullfile(basepath, dirnames(i)));
        if ~exist(filepath, 'dir')
            warning('%s does not exist, skipping...', filepath)
            continue
        end
        cd(filepath)
        
        for ii = 1 : length(vars)           
            filename = dir(['*', vars{ii}, '*']);
            if length(filename) == 1
                d{i, ii} = load(filename.name);
            else
                warning('no %s file in %s, skipping', vars{ii}, filepath)
            end
        end
    end
end

% params
% spkgrp = f{1}.spkgrp;
% ngrp = length(spkgrp);
% nsessions = length(f);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% analyze data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fr = cell(1, nsessions);
if forceA
    for i = 1 : nsessions
        filepath = char(fullfile(basepath, dirnames(i)));
        cd(filepath)
        [~, basename] = fileparts(filepath);
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % session info (cell explorer foramt)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        session = sessionTemplate(pwd, 'showGUI', false);
        cell_metrics = ProcessCellMetrics('session', session,...
            'manualAdjustMonoSyn', false, 'summaryFigures', false,...
            'debugMode', true, 'transferFilesFromClusterpath', false,...
            'submitToDatabase', false);
        
        
        if ~isempty(d{i, 2})
            spikes = d{i, 2}.spikes;
            fr{i} = FR(spikes.times, 'basepath', filepath, 'graphics', false, 'saveFig', false,...
                'binsize', 60, 'saveVar', false, 'smet', 'MA', 'winBL', [1 Inf]);
        end
        
        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rearrange data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all

% fr
figure
for i = 1 : nsessions
    subplot(1, nsessions, i)
    if ~isempty(fr{i})
stdshade(fr{i}.strd(:, :), 0.3, 'k', fr{i}.tstamps / 60, 3)
ylim([0 20])
        %         plot(fr{i}.norm')
    end
end

% states
load([basename '.SleepState.states.mat'])
wake = SleepState.ints.WAKEstate;
figure
plot(acc.tstamps, acc.data)
hold on
axis tight
y = ylim;
plot([wake(:, 1) wake(:, 1)], y, 'k')
plot([wake(:, 2) wake(:, 2)], y, 'r')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% to prism
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
