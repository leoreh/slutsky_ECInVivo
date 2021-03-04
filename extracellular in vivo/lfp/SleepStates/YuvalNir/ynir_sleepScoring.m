function sleepScoring (EEG1, EEG2, EMG, existing_scoring, right_scoring)
close all;  tic; maxLen = length(EMG); 
if ~exist('right_scoring','var'), right_scoring = []; end 
%% 0) User Parameters
    SAMPLE_RATE = 1000; % WARNING: NOT working in another rate other then 1000 Hz!   <------------
    EEG_STD_4_YLIM = 15;
    EMG_STD_4_YLIM = 7.5;
    
    SCORING_FILE_NAME = 'sleep_scoring'; % Where to save files? absoulut path or file name in current directory

    Active_Window_Size = 10; % (seconds) Active window size - Scoring was done in 4 second epochs
    Inactive_Wing_Size = 2; % (seconds) size of pre & post window activity    
	automaticSaveEvery_X_minutes = nan; % Make nan or inf to disable
    
% Keybord hot keys
%--------------------------------------------------------------------------
    ENTER = 13;     % Enter
    LEFT_KEY = 52;  % The number 4
    RIGHT_KEY = 54;	% The number 6
    I_KEY = 105;    % The letter 'i' (For information)
    Z_KEY = 122;    % The letter 'z' (Change hipnogram - relative/absolute)
    T_KEY = 116;    % The letter 't' (Show "master" scoring)
    C_KEY = 99;     % The letter 'c' (Check scoring)
    A_KEY = 97;     % The letter 'a' (Automated score into master scoring)
    PLUS_KEY = 43;  % The operation '+'
    MINUS_KEY = 45;	% The operation '-'
    SHIFT_PLUS_KEY = 91;  % The sign '['
    SHIFT_MINUS_KEY = 93;	% The sign ']'
                    
    % Scoring Keys =    [   0          1           2         3           5        7        8      9         '/'            ' *'          
    %                   -------------------------------------------------------------------------------------------------------------           
    scoreKeys =         [  48   ,     49,         50,       51,         53,      55,      56,    57,        47,             42       ]; 
    scroeStrings =      {'Event','ActiveWake',  'ARTIFACT','ARTIFACT3', 'UNKNOWN', 'WAKE', 'NREM', 'REM','TRANSITION-WN','TRANSITION-NR'};
    outputScoreValues = [   0,       -1,           -2,        -3,           2,       4,       5,     6,        9,            10       ];
    % Test char keys: ### isAKey = waitforbuttonpress; double(get(gcf,'CurrentCharacter')), close(gcf); ### 

    % Auto complete states
    autoCompleteStates = {'Event', 'ARTIFACT'}; % Change only this line, not the one below (which automatically extract the corresponding values).


%% 1) Automatic pre-processing and initializations  

% Self Initializtions - not for user
cX = []; 
autoCompleteStartTime = 0; 
hypnoZoomFlag =1; showTeacherFlag = 1;  
AUTO_COMPLETE_VALUES = nan(size(autoCompleteStates)); for autoCompleteI = 1:length(autoCompleteStates), AUTO_COMPLETE_VALUES(autoCompleteI) = score_str2value(autoCompleteStates(autoCompleteI)); end
if ~exist('existing_scoring','var') || isempty(existing_scoring), existing_scoring = ones(1,length(EEG1))*score_str2value('UNKNOWN'); end
%-----------------------------------

[problematicEEG,EEG1, EEG2] = checkAndRestoreSignals(EEG1, EEG2); % Fix marker of noisy EEG channel (NaN at start of vector)

WINDOW_SIZE = Active_Window_Size * SAMPLE_RATE; %(ms) 
WING_SIZE = Inactive_Wing_Size * SAMPLE_RATE;   %(ms) 
i = WINDOW_SIZE + WING_SIZE;                        % The prime variable for iteration number
fig = figure;  fullScreen(fig);  hZoom = zoom(fig); % Open figure, enlarge & define zoom mode
REMOVAL_RADIUS = WINDOW_SIZE /66.667; %(ms) How far (in ms) from a mark you should press, to remove it
sigLen = length(EEG1);
axisMemory = [];

% Y-Limits of EEG & EMG plots
    eeg_avg = mean([EEG1,EEG2]); eeg_std = std([EEG1,EEG2]);
    emg_avg = mean(EMG); emg_std = std(EMG);
    yLimEEG = changeYlim(eeg_avg, eeg_std, EEG_STD_4_YLIM, 1);
    yLimEMG = changeYlim(emg_avg, emg_std, EMG_STD_4_YLIM, 1); 

% Scoring lists, & Scoring data : If available - translate, else - vectors of nans  
    [scoring, scoreTiming] = explicitScoring(existing_scoring, 1);              % Exisiting
    if ~isempty(right_scoring) 
        [right_scoring, right_scoreTiming] = explicitScoring(right_scoring);    % Right
    else, right_scoreTiming = []; end
    

    
%% 2) Run Sleep Scoring Algorithm
try
while true
    automaticSave(automaticSaveEvery_X_minutes, scoring, scoreTiming);
    
    %% A) Draw Graphs
    x = (i-(WINDOW_SIZE + WING_SIZE)+1) : (i+WING_SIZE);   % Absolute time
    relativeX = (1:length(x))-WING_SIZE;                            % Relative time
    eeg1 = EEG1(x);     eeg2 = EEG2(x);
    
    % 1) Plot Hypnogram
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    gca1 = subplot(24,1,1:2);   cla;
    
    % Plot teacher hypnogram (If availabel)
    if ~isempty(right_scoreTiming) && showTeacherFlag
        [timesF, scoreF, timeEvents, scoreEvents] = scoringToHypongram(right_scoreTiming, right_scoring);
        if ~isempty(timeEvents), hold on; plot(timeEvents, scoreEvents,'xr'); hold off; end; % Plot artifacts and events
        line(timesF, scoreF, 'color', [1,0,0]);                                 % Plot graph        
    end
    
    % Plot hypnogram 
    [timesF, scoreF, timeEvents, scoreEvents, lastScoreTime] = scoringToHypongram(scoreTiming, scoring);
    if lastScoreTime < x(end), timesF = [timesF, x(end)]; scoreF = [scoreF, scoreF(end)]; end % If we are in a new "area" (no markers further (yet)) - insert too more coordinates so the current stage will appear in the hypnogram "online"
    if ~isempty(timeEvents), hold on; plot(timeEvents, scoreEvents,'xc'); hold off; end; % Plot artifacts and events
    line(timesF, scoreF, 'color', [0,0,0]);                                 % Plot graph
    if hypnoZoomFlag, hypnoXlim = max(timesF(end),lastScoreTime); else hypnoXlim = length(EEG1); end; 
    
    % Design
    line([1,1]*x(1)+WING_SIZE, [score_str2value('UNKNOWN')-1, score_str2value('REM')+1], 'color',[1,0,0]);           % Plot current point line
    set(gca,'xtick',[],'yTick',[score_str2value('UNKNOWN'), score_str2value('WAKE'), score_str2value('NREM'), score_str2value('REM')])
    set(gca,'yTickLabel',{'UNKNOWN','WAKE','NREM','REM'},'TickLength',[0,0.01],'fontsize',8);   
    xlim([0, hypnoXlim]); ylim([score_str2value('UNKNOWN')-1, score_str2value('REM')+1]); line([0, hypnoXlim], [1,1]*54, 'color', [1,1,1]*0.9);
    tit = sprintf('Progress: %.2f%% (%s hours).    For more information, options and shortcuts press  ''i'' ',(100*i/sigLen), datestr((i/1000-10)/(24*60*60), 'HH:MM:SS')); title(tit,'fontsize',11);
    
    
    % 2) Plot Signal Graph
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    gca2 = subplot(24,1,3:14); stLnCl = 0.2;
    plot (relativeX([1,end]),[1,1]*(+1)*diff(yLimEEG)/6,'color',[1,1,1]*stLnCl); line(relativeX, eeg1+diff(yLimEEG)/6,'color',[0,0,0.7]); 
    line (relativeX([1,end]),[1,1]*(-1)*diff(yLimEEG)/6,'color',[1,1,1]*stLnCl); line(relativeX, eeg2-diff(yLimEEG)/6,'color',[0.2,0.4,0]); 
    ylim(yLimEEG); xlim([relativeX(1), relativeX(end)]);  ylabel('Voltage (uV)'); % Plot relative time 
    set(gca2,'xtick',((relativeX(1):1000:relativeX(end))-1)); set(gca2,'xticklabel',((relativeX(1):1000:relativeX(end))-1)/1000);
    
    % Plot already marked scoring
    currentScroingsIdx = (find(x(1) < scoreTiming & scoreTiming < x(end)));
    for si = 1:length(currentScroingsIdx)
          markScore(scoreTiming(currentScroingsIdx(si))-x(1), scoring(currentScroingsIdx(si))); % plot marker and text on graph
    end
    if x(1)>1, markScore(0, scoring(scoreTiming == max(scoreTiming(scoreTiming < x(1)))), 1); end        % If there is a prev window, Plot previous mark

    % Plot Right Scoring
    if ~isempty(right_scoreTiming) && showTeacherFlag
        currentScroingsIdx = (find(x(1) < right_scoreTiming & right_scoreTiming < x(end)));
        for si = 1:length(currentScroingsIdx)
            markScore(right_scoreTiming(currentScroingsIdx(si))-x(1), right_scoring(currentScroingsIdx(si)), 0, 1); % plot marker and text on graph
        end
        if x(1)>1,  markScore(0, right_scoring(right_scoreTiming == max(right_scoreTiming(right_scoreTiming < x(1)))), 1, 1); end  % If there is a prev window, Plot previous mark
    end       

    % Gray wing areas
    plotBackground([0,WING_SIZE-1]+relativeX(1),0,yLimEEG(2),'k'); plotBackground([-WING_SIZE+1,0]+relativeX(end),0,yLimEEG(2),'k');                 % Plot the gray background

    
    % 3) Plot EMG
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    gca3 = subplot(24,1,15:19);
    plot(relativeX, EMG(x),'color',[1,0.5,0]);   ylim(yLimEMG); xlim([relativeX(1), relativeX(end)]);
    plotBackground([0,WING_SIZE-1]+relativeX(1),0,yLimEMG(2),'k'); plotBackground([-WING_SIZE+1,0]+relativeX(end),0,yLimEMG(2),'k');                 % Plot the gray background
    set(gca3,'xtick',((relativeX(1):1000:relativeX(end))-1)); set(gca3,'xticklabel',((relativeX(1):1000:relativeX(end))-1)/1000);

    
    % 4) Plot Spectrogram
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    gca4 = subplot(24,1,20:24);    
    [~,F,T,P] = spectrogram(eeg1, 750, 0, 0:0.3:20, SAMPLE_RATE, 'yaxis'); %colormap(jet);
%     [~,F,T,P] = spectrogram(eeg1, 2000, 1000, 0:0.3:20, SAMPLE_RATE, 'yaxis'); colormap(gray); % Amit's idea
    P= P./repmat(max(P,[],1),[size(P,1), 1]); % Normalize columns
    T = 1000*T - WING_SIZE; imagesc( T, F, (P) );  set(gca,'YDir', 'normal'); 
    ylabel('Normalized Hz'); xlabel('Time (s)');  xlim([relativeX(1),relativeX(end)]);
    line(get(gca,'xlim'), [1,1]*5, 'color' ,[1,1,1]);
    ylim([0,20]); caxis([min(P(:)), max(P(:))*1.15]);
    
    set(gca4,'xtick',(( round(relativeX(1)/1000)*1000:1000:relativeX(end)))); set(gca4,'xticklabel',(( round(relativeX(1)/1000)*1000:1000:relativeX(end)))/1000);
    hold off
    
    %% B) User Interaction
    while true 
        isAKey = waitforbuttonpress;                                        % Wait for user (mouse / keyboard)
        key = double(get(fig,'CurrentCharacter'));                          % Last key pressed (even if moused was pressed in the previous line)
        if gcf ~= fig, continue; end;                                       % If input not coming from main window - continue

        % 1) If mouse press
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ~isAKey()                                            
            cor = double(get(gca,'CurrentPoint')); cX = cor(1)+x(1)+WING_SIZE;      % Get X coordinate of the press
            axisMemory = gca;
            
            % Validity checks for mouse press
            if (gca == gca1) && strcmp(hZoom.Enable,'off'),                 % If pressed in hypnogram (& not in zoom mode), 
                newI = WING_SIZE + ceil(cor(1)/WINDOW_SIZE)*WINDOW_SIZE;    % jump to place, and break to plot all graphs
                if 0 < newI - WING_SIZE && newI + WING_SIZE <= sigLen, i = newI; break; end;  
            elseif (gca == gca2) && ( cor(1)<0 || WINDOW_SIZE<cor(1) )      % if press in EEG graph was indside a wing ... 
                cX = []; continue;                                          %       Cancel press and get input again
            end
            
            % Marker manipulation with only mouse-
            % If pressed near previous marks, delete them - delete from lists and break to plot again
            delSpIdx = find(((scoreTiming-cX)/REMOVAL_RADIUS).^2 <= 1); % Find markers in 'REMOVAL_RADIUS' from the X coordinate
            if ~isempty(delSpIdx) && (gca == gca2), scoreTiming(delSpIdx)= nan; scoring(delSpIdx)= nan;  break;
            elseif autoCompleteStartTime > 0 && (gca == gca2)       % Auto-complete last score
                previousMarkIdx = max(find(scoreTiming < autoCompleteStartTime)); 
                hold on; markScore(cX-x(1), scoring(previousMarkIdx)); hold off;                         % plot marker and text on graph
                [scoring, scoreTiming] = implementScore(scoreTiming, scoring, cX, scoring(previousMarkIdx));    % Add score to the 2 lists
                autoCompleteStartTime = 0;  break;   % Auto was completed- init. 
            end
            
        % 2) If keyboard pressed
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        else
            if ~isempty(axisMemory), axes(axisMemory); end                  % If want to do the same action more then once (e.g. '+'), then after first time, the code will break to draw all over again, and current axis will be lost - so load it [= current axis chages only by mouse bottom].
            switch key
            case ENTER,         hereSave(scoring, scoreTiming,1);  % Save lists of scorings
            case I_KEY,         instructionWindow(); set(0, 'currentfigure', fig);  
            case Z_KEY,         hypnoZoomFlag   = 1 - hypnoZoomFlag;   
            case T_KEY,         showTeacherFlag = 1 - showTeacherFlag;  
            case C_KEY,         checkScoringHere(EEG1, EEG2, EMG, scoring, scoreTiming); set(0, 'currentfigure', fig);  
            case A_KEY,         [right_scoring, right_scoreTiming] = runAutomatedScoring(EEG1, EEG2, EMG, SAMPLE_RATE);
            case RIGHT_KEY,   i = nextStep(i,+1);   
            case LEFT_KEY,   i = nextStep(i,-1);  

            case PLUS_KEY 
                if gca == gca3, [yLimEMG, EMG_STD_4_YLIM] = changeYlim(emg_avg, emg_std, EMG_STD_4_YLIM, 0.5); 
                else [yLimEEG, EEG_STD_4_YLIM] = changeYlim(eeg_avg, eeg_std, EEG_STD_4_YLIM, 0.5);  end
            case MINUS_KEY
                if gca == gca3, [yLimEMG, EMG_STD_4_YLIM] = changeYlim(emg_avg, emg_std, EMG_STD_4_YLIM, 2); 
                else [yLimEEG, EEG_STD_4_YLIM] = changeYlim(eeg_avg, eeg_std, EEG_STD_4_YLIM, 2);     end  
                                
            case SHIFT_PLUS_KEY 
                [i, WINDOW_SIZE, WING_SIZE] = changeXaxis(-1, i, WINDOW_SIZE, WING_SIZE);
                
            case SHIFT_MINUS_KEY
                [i, WINDOW_SIZE, WING_SIZE] = changeXaxis(+1, i, WINDOW_SIZE, WING_SIZE);
                
            otherwise
                if ismember(key,scoreKeys) && (gca == gca2) && ~isempty(cX)
                    scoreV = score_key2Value(key);
                    if ismember(scoreV, AUTO_COMPLETE_VALUES), autoCompleteStartTime = cX; end    % If key need an Auto complete, keep location press
                    hold on; markScore(cX-x(1), scoreV); hold off;                             % plot marker and text on graph
                    [scoring, scoreTiming] = implementScore(scoreTiming, scoring, cX, scoreV); % Add score to the 2 lists
                else
                    continue;   % If the key is not a possible key - start input loop again
                end
            end
            break; % After resolving key action - break input loop, and plot all again  % break the inner while loop (User Interaction)  % After each mark, plot everything again (bottom line - to update hypnogram)
        end
    end  % Only keyboard press can break the while
    cX = [];
end
catch ME,   
    if strcmp(ME.message, 'waitforbuttonpress exit because target figure has been deleted')
         hereSave(scoring, scoreTiming, 1);
    else
        hereSave(scoring, scoreTiming, 1); disp('@@@@@@@@@@@@@@@@@@@@@@@@@@@@');
        disp('!! PROGRAM CRASHED!! scoring was automatically saved!'); disp('Error is:'); throw(ME); 
    end    
end

%% 3) ~~~~~~~~~ Inner Functions ~~~~~~~~~ 

function [problematicEEG,EEG1, EEG2] = checkAndRestoreSignals(EEG1, EEG2)
    if isnan(EEG1(1))
        EEG1(1) = 0; problematicEEG = 1;
    elseif isnan(EEG2(1))
        EEG2(1) = 0; problematicEEG = 2;
    else
        problematicEEG = 0;
    end
end

function [lim, std4lim] = changeYlim(avg, std, std4lim, change)
    std4lim = std4lim * change; lim = [-1,1] * (avg + std * std4lim);        
end

function [timesF, scoreF, timeEvents, scoreEvents, lastScoreTime] = scoringToHypongram(scoreTiming, scoring)
    chosenIdxs = (~isnan(scoreTiming));
    times = scoreTiming(chosenIdxs);      scores = scoring(chosenIdxs);      % Remove Nans
    [times , idxs] = sort(times);   scores = scores(idxs);                    % Sort
    chosenIdxs = ismember(scores, AUTO_COMPLETE_VALUES);
%     chosenIdxs = ( scores == score_str2value('Event') ) | ( scores == score_str2value('ARTIFACT') );
    timesF = times(~chosenIdxs);              scoreF = scores(~chosenIdxs);  % Find artifacts and events

    len = length(timesF); % Duplicate points for straight lines, and add nans to remove vertical lines
    timesF = reshape([[nan,timesF(2:end)]; nan(1,len); timesF], 1, 3*len);
    scoreF = reshape([nan(1,len); scoreF; [scoreF(1:end-1),nan]], 1, 3*len);
    timesF = timesF(3:end);     scoreF = scoreF(2:end-1); 

    timeEvents = times(chosenIdxs); scoreEvents = scores(find(chosenIdxs)-1);

    lastScoreTime = times(end);
end
                
function [i, WINDOW_SIZE, WING_SIZE] = changeXaxis(sign, i, WINDOW_SIZE, WING_SIZE)
    dWind = sign*500; dWing = round((WINDOW_SIZE + sign*500)*WING_SIZE/WINDOW_SIZE) - WING_SIZE;

    if 1000 - WINDOW_SIZE <= dWind && dWind <= sigLen - i - WING_SIZE - dWing % [WINDOW_SIZE+dWind >= 1000] To block zoom-in more then 1 sec. [i + dWind + WING_SIZE + dWing <= sigLen] To block length more then signal len
        WINDOW_SIZE = WINDOW_SIZE + dWind; WING_SIZE = WING_SIZE+dWing; i = i + dWind + dWing;
    end
end

function markScore(x_of_fig, score, pastMarkFlag, instructiveFlag)
    if nargin < 3, pastMarkFlag = 0; end
    if nargin < 4, instructiveFlag = 0; end
    
    color = [0,1,0]*instructiveFlag + [1,0,0]*(1-instructiveFlag);
    
    if ~pastMarkFlag, line([1,1]*(x_of_fig-WING_SIZE), yLimEEG, 'color', color, 'linewidth',2); end % Plot Line
    
    if instructiveFlag, y = yLimEEG(2)-diff(yLimEEG)/100; alignArg = 'right';
    else                y = yLimEEG(1)+diff(yLimEEG)/100; alignArg = 'left'; end
    
    tempCord = double(round([(x_of_fig-WING_SIZE)+(WING_SIZE*2+WINDOW_SIZE)/200, y]));
    text(tempCord(1), tempCord(2), score_value2str(score) , 'color', [0.5,0,1], 'FontSize',13, ...
        'HorizontalAlignment', alignArg, 'FontWeight', 'bold', 'Parent', gca2,'Rotation',90);
end

function [scoring, scoreTiming] = implementScore(scoreTiming, scoring, time, score) % Add score to the 2 lists
    tmp_idx = min(find(isnan(scoreTiming)));  scoreTiming(tmp_idx) = time;    scoring(tmp_idx) = score; 
end

function hereSave(scoring, scoreTiming, timeFlag), if nargin < 3, timeFlag = 0; end
    [~, tmp_idx] = sort(scoreTiming); scoreTiming = scoreTiming(tmp_idx); scoring = scoring(tmp_idx); % Sort; 
    
    scoring = compactScoring(scoring, scoreTiming);
    scoring = scoring(1:min(maxLen,length(scoring)));

    if timeFlag, save([SCORING_FILE_NAME,'_',datestr(now,30)], 'scoring'); 
    else save([SCORING_FILE_NAME,'_AutoSave'], 'scoring'); end
end

function i = nextStep(i, direction)
    
    new_i = i + direction * (WINDOW_SIZE);   
    
    if new_i < WINDOW_SIZE + WING_SIZE || sigLen < new_i + WING_SIZE, return; end;
    
    i = new_i;        
end

function instructionWindow()
   iFig = figure;
   set(gca,'Visible','off')
   fS = 14;   
   ls = 0.12; %lineSpace
   
      
   controlRIGHT = '"6": move forward.';
   controlLEFT = '"4": move backward.';
   zoomX = '"+" or "-": zooming in/out in the y (amplitude) axis';
   zoomY = '"[" or "]": zooming in/out in the x (time) axis';
   hypno = '"z": change relative/absolute hypnogram'; 
   guide = '"t": reveal/hide guiding socring (if available)';
   check = ' "c": check scoring by scatter distribution';
   enter = '"Enter": saves current progress in current directory';
   ScoringIns = 'To score press left click & choose a scoring key';
   ScoringKeys = '(scoring keys are 0,1,2,3,5,7,8,9,/,*).';
   

    text(0.35, 1, 'Instructions','Units','normalized', 'FontWeight', 'bold', 'HorizontalAlignment','left', 'FontSize',fS);
    text(0, 0.85, controlRIGHT,'Units','normalized', 'HorizontalAlignment','left', 'FontSize',fS);
    text(0, 0.85-0.5*ls, controlLEFT,'Units','normalized', 'HorizontalAlignment','left', 'FontSize',fS);
    text(0, 0.85-1.5*ls, zoomX,'Units','normalized', 'HorizontalAlignment','left', 'FontSize',fS);
    text(0, 0.85-2*ls, zoomY,'Units','normalized', 'HorizontalAlignment','left', 'FontSize',fS);
    text(0, 0.85-3*ls, hypno,'Units','normalized', 'HorizontalAlignment','left', 'FontSize',fS);
    text(0, 0.85-4*ls, guide,'Units','normalized', 'HorizontalAlignment','left', 'FontSize',fS);
    text(0, 0.85-5*ls, check,'Units','normalized', 'HorizontalAlignment','left', 'FontSize',fS);
    text(0, 0.85-6*ls, enter,'Units','normalized', 'HorizontalAlignment','left', 'FontSize',fS);
    text(0, 0.85-7*ls, ScoringIns,'Units','normalized', 'HorizontalAlignment','left', 'FontSize',fS);
    text(0, 0.85-7.5*ls, ScoringKeys,'Units','normalized', 'HorizontalAlignment','left', 'FontSize',fS);
end

function v = score_key2Value(k), v = outputScoreValues(k == scoreKeys); end
function s = score_value2str(v), s = scroeStrings{v == outputScoreValues}; end
function v = score_str2value(s), v = outputScoreValues(strcmp(scroeStrings,s)); end

%% 3) ~~~~~~~~~ Imported Functions ~~~~~~~~~ 

function plotBackground( X, Y, DEV, c)
    % Plot Backround around graph Y=f(X) in size of (vector) DEV and color c.
    % If Y is a single value (length=1) the backround will not follow the
    % graph. For example put Y = 0 to make the backround mark the specific area
    % on the graph.
    % Default color is blue ('b').

    if nargin < 4
       c = 'b'; 
    end

    if length(Y) == 1
        Y = ones(size(X))*Y;
    end

    ebars=patch([X,fliplr(X)],[Y-DEV,fliplr(Y+DEV)],c);
    set(ebars,'EdgeColor','none'); % Change edges color.

    alpha(0.25);
end



function fullScreen(f), set(f,'units','normalized','outerposition',[0 0 1 1]); end

function checkScoringHere(eeg1, eeg2, emg, scoringCMP, scoringTiming)        
    sleepScoringTesting(eeg1, eeg2, emg, compactScoring(scoringCMP, scoringTiming))    
end

    
    
function scoring = compactScoring(scoring, scoreTiming)
    lastScoreIdx = max(find(~isnan(scoreTiming)));   scoreTiming = round(scoreTiming);    % Find last index of scoring & round timing to integers (ms)
    scoring(lastScoreIdx) = score_str2value('UNKNOWN');                             % Set last score to UNKNOWN
    scoring(scoreTiming(lastScoreIdx):end) = scoring(lastScoreIdx);                 % Fill all until end with new "last score" - UNKNOWN
    for idx = (lastScoreIdx-1):-1:1                                                 % Fill (backward) all array with scoring in ms resolution
        scoring(scoreTiming(idx):(scoreTiming(idx+1)-1)) = scoring(idx); 
    end
end

function [scoring, scoreTiming] = explicitScoring(scoring, nanPadFlag)  % True only for 1000 Hz
    if nargin < 2, nanPadFlag = 0; end
    
    if iscolumn(scoring), scoring = scoring'; end
    scoringChangeIndexes = [1,find(diff(scoring) ~= 0)+1];
    
    if nanPadFlag
        nanPadding = nan(1,length(scoring) - length(scoringChangeIndexes)+1);
    else
        nanPadding = [];
    end
    
    scoreTiming = [    scoringChangeIndexes  , nanPadding];
    scoring = [scoring(scoringChangeIndexes) , nanPadding];
    
	if isnan(scoring(1)),  scoring(1) = score_str2value('UNKNOWN'); scoreTiming(1) = 1; end % Always put UNKNOWN at start of NEW scoring                
end

function [auto_scoring, auto_scoreTiming] = runAutomatedScoring(EEG1, EEG2, EMG, SAMPLE_RATE)
    % Call automatic scoring and present it insted of teacher (green) scoring
    
    % Currently not working
    auto_scoring = []; auto_scoreTiming = [];    return
    
    tic
    try
        fullScoring = nirScoring([EEG1; EEG2; EMG], SAMPLE_RATE, 0.8);
        [auto_scoring, auto_scoreTiming] = explicitScoring(fullScoring);
    catch
        auto_scoring = []; auto_scoreTiming = []; return;
    end
    toc
end

function automaticSave(automaticSaveEvery_X_minutes, scoring, scoreTiming)
   if toc >= automaticSaveEvery_X_minutes * 60
       f = waitbar(0.5,'Automatic saving...');
       hereSave(scoring, scoreTiming, 0);
       close(f);
       tic
   end
end

end