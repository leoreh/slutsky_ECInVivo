% states_sessions

% organizes and plots states

forceA = false;
forceL = false;
saveFig = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

basepath = 'F:\Data\Processed\lh58';
vars = ["session.mat";...
    "SleepState.states.mat"];      
pcond = ["tempFlag"];     
ncond = [""];                      

if ~exist('varArray', 'var') && ~forceL
    [varArray, dirnames] = getSessionVars('basepath', basepath, 'vars', vars,...
        'pcond', pcond, 'ncond', ncond, 'sortDir', true);
end
nsessions = length(dirnames);

session = varArray{1, 1}.session;
ss = varArray{1, 2}.SleepState;
nchans = session.extracellular.nChannels;
fs = session.extracellular.sr;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% analyze data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if forceA
    close all
    for i = 1 : nsessions
       
        % file
        filepath = char(fullfile(basepath, dirnames(i)));
        cd(filepath)
        
        % session info
        session = CE_sessionTemplate(pwd, 'viaGUI', false,...
            'force', true, 'saveVar', true);
        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rearrange data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1 : nsessions
    
    % file
    ss = varArray{i, 2}.SleepState;
    
    states = {ss.ints.WAKEstate, ss.ints.NREMstate, ss.ints.REMstate};
    recDur(i) = max(max(vertcat(states{:})));
    for ii = 1 : length(states)
        stateDur(i, ii) = sum(states{ii}(:, 2) - states{ii}(:, 1))...
            / recDur(i) * 100;
    end
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all

pathPieces = regexp(dirnames(:), '_', 'split'); % assumes filename structure: animal_date_time
sessionDate = [pathPieces{:}];
sessionDate = sessionDate(2 : 3 : end);

pflag = 1;
grp = [1 : 8];
if pflag
    fh = figure;
    bar(stateDur, 'stacked')
    legend({"Wake"; "NREM"; "REM"})
    xticks(1 : nsessions)
    xticklabels(sessionDate)
    title('Sleep States')
    xlabel('Session')
    ylabel('Relative Duration [%]')
    box off
    if saveFig
        figname = fullfile(basepath, 'UnitsDetected');
        % print(fh, figname, '-dpdf', '-bestfit', '-painters');
        export_fig(figname, '-tif', '-transparent', '-r300')
    end
end