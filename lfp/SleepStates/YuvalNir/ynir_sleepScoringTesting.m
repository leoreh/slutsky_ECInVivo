function sleepScoringTesting(eeg1, eeg2, emg, scoring)
%% Define
fig = figure; fullScreen(fig); 

freqs = [0:0.2:100]; % Relevant freqencies for power analysis
normalizedSummaryScatterPlot = 0;


% Constant variables
states = {'Active', 'Wake', 'N1', 'NREM', 'REM'}; 
colors = [[1,1,0]*0.9; 1,0,0; 0,0.75,0.75; 0,0,1; 0,1,0]; 
statesCode = [-1, 4,9,5,6];
hypnoText = {'UNKWONW','Wake', 'NREM', 'REM'};  
hypnoStateCode = [2,4,5,6];  
hypnoTwinState =[-1, 4; 9, 5; 10, 6]; % Leave empty "[]" for ignoring twin states

subplots = {1:12, [13:15,25:27], [37:39,49:51], [61:63,73:75,85:87],[78,79,90,91], [18,19,30,31,42,43,54,55,66,67], [80:84,92:96], [20:24,32:36,44:48,56:60,68:72]};
plotHypno=1; plotFFT1=2; plotFFT2=3; plotBars=4; plotLegend=5; plotEMG=6; plotEEG=7; plotScatter=8;
plotRows = 8; plotCols=12;
nStates = length(states);

scoring = scoring(1:(max(find(diff([0,scoring])))+1)); % Cut all padding 'UNKNOWN's at the end
fWB = waitbar(0,'Analyzing Scoring...'); percents = 0; percentChange = 1/nStates/4;
%% Main
for state_i = 1:nStates
    c = colors(state_i, :);
    timeRanges = splitedTimings (scoring, statesCode(state_i), 4000, 2000); % State timing
    waitbar(percents, fWB, 'Analyzing Scoring...'); percents = percents + percentChange;
    
    [outputs1, fftFP1] = powerRmsAndSpectrum(eeg1, emg, timeRanges);    % EEG, EMG
    waitbar(percents, fWB, 'Analyzing Scoring...'); percents = percents + percentChange;
    [outputs2, fftFP2] = powerRmsAndSpectrum(eeg2, emg, timeRanges);
    waitbar(percents, fWB, 'Analyzing Scoring...'); percents = percents + percentChange;
    plotStateEEG_EMG(outputs1, fftFP1, outputs2, fftFP2, c);
    waitbar(percents, fWB, 'Analyzing Scoring...'); percents = percents + percentChange;
end
close (fWB);

legendPlot(states, colors);
plotPercentageState(scoring, states, statesCode, colors);
scoringToHypongram(scoring, hypnoStateCode, hypnoText, hypnoTwinState)
designSubplots();

% SAVE
% elabel= ['scoring_E', num2str(E)];
% saveas(gcf, [pathTo_SDC('scoringOutputFolder', E), '\',elabel, '.jpg']);    % Save JPG
% saveas(gcf, [pathTo_SDC('scoringOutputFolder', E), '\',elabel]);            % Save Figure


%% Functions
% ~~~~~~~~~~~~~~~~~~~
function timeRanges = splitedTimings (scoring, state, size, minLimit)
    stateScoring = single(scoring == state);
    
    for i = 2:length(stateScoring)
        stateScoring(i) = (stateScoring(i-1)+stateScoring(i)) * stateScoring(i);
    end
    
    cutStateIdx = (mod(stateScoring+1,size+1) == 0);    % find where to cut
    stateScoring = scoring == state;        % get binary vector again
    stateScoring(cutStateIdx) = 0;          % zero the cut indexes
    
    d = diff([0, stateScoring, 0])';                 % Translate scoring
    timeRanges = [find(d == 1), find(d == -1)-1];       % Find start & end
    timeRanges (diff(timeRanges,[],2) < minLimit,:) = [];   % Remove too short
end

function [outputs, fftFP] = powerRmsAndSpectrum(eeg, emg, stateTiming)

    if isnan(eeg(1)), fftFP = nan(2,1); outputs = []; return; end;
    if isempty(stateTiming), fftFP = nan(2,1); outputs = []; return; end;
    
    outputs = nan(length(stateTiming), 2); sumfft = 0;
    for ii = 1:length(stateTiming)
        % Range
        range = stateTiming(ii,1):stateTiming(ii,2);
        % FFT
        [pow,f] = pwelch(eeg(range),round(1000)*2,round(1000*0.25),freqs, 1000);
        % Save power & rms
        outputs(ii,1) = sum(pow(f > 25)) / sum(pow(f < 5));
        outputs(ii,2) = sqrt(mean(emg(range).^2));
        sumfft = sumfft + pow; 
    end 
    outputs = log(outputs);
    fftFP = [f; sumfft/length(stateTiming)];
end

function plotStateEEG_EMG(outputs1, fftFP1, outputs2, fftFP2, stateColor)

    subplot(plotRows,plotCols,subplots{plotFFT1}); hold on    % plotStateSpectrum1
    plot(fftFP1(1,:) ,log(fftFP1(2,:)), 'color', stateColor); 
    
    subplot(plotRows,plotCols,subplots{plotFFT2}); hold on    % plotStateSpectrum2
    plot(fftFP2(1,:) ,log(fftFP2(2,:)), 'color', stateColor); 
    
    outputs = [outputs1 ; outputs2]; if isempty(outputs), return; end
    subplot(plotRows,plotCols,subplots{plotScatter});  hold on   % plotMainScatter
    plot(outputs(:,1), outputs(:,2), '.', 'color', stateColor);

    normFlag = normalizedSummaryScatterPlot;    
    subplot(plotRows,plotCols,subplots{plotEEG}); hold on    % plotSubRatio
    [tmpb,tmpa] = histcounts(outputs(:,1),'BinWidth',0.1) ;
    plot((tmpa(1:end-1)+tmpa(2:end))/2, tmpb/(max(tmpb)).^normFlag, 'color', stateColor, 'linewidth',2);
    
    subplot(plotRows,plotCols,subplots{plotEMG}); hold on 	% plotSubRMS
    [tmpb,tmpa] = histcounts(outputs(:,2),'BinWidth',0.1) ;
    plot((tmpa(1:end-1)+tmpa(2:end))/2, tmpb/(max(tmpb)).^normFlag, 'color', stateColor, 'linewidth',2);    
end

function legendPlot(text, colors)
    
    text = {'~~ LEGEND ~~', text{:}};
    colors = [1,1,1; colors];
    
    subplot(plotRows,plotCols,subplots{plotLegend}); 
    for i = 1:length(text)
        hold on; plot([0,0],[0,0], 'color', colors(i,:),'linewidth',3); hold off;
    end
    legend(text{:},'Location','southwest');
    axis off;
end

% scoringToHypongram(scoring, [2,4,5,6], {'UNKNOWN','WAKE','NREM','REM'} , [100,-1,9,10])
function scoringToHypongram(scoring, mainScores, titles, semiScores)
    % semiScores will appear as mainScores but softer color (e.g. transitions)
    
    subplot(plotRows,plotCols,subplots{plotHypno});
    d = diff([0, scoring]);                 % Translate scoring
    timesF = find(d ~= 0);
    scoreF = scoring(timesF);

    len = length(timesF); % Duplicate points for straight lines, and add nans to remove vertical lines
    timesF = reshape([[nan,timesF(2:end)]; nan(1,len); timesF], 1, 3*len);
    scoreF = reshape([nan(1,len); scoreF; [scoreF(1:end-1),nan]], 1, 3*len);
    timesF = timesF(3:end);     scoreF = scoreF(2:end-1); 

    
    tmp = zeros(size(scoreF));
    for i = 1:length(mainScores)
        tmp = tmp + (scoreF == mainScores(i));
    end
    for i = 1:size(semiScores,1)
        tmp = tmp + (scoreF == semiScores(i,1))*2;
    end
    mainIdx = (tmp == 1); semiIdx = (tmp == 2);
    
    mainS = scoreF; mainS(semiIdx) = nan;
    semiS = scoreF; semiS(mainIdx) = nan;
    
    
    for i = 1:size(semiScores,1)
        semiS(semiS == semiScores(i,1)) = semiScores(i,2);
    end
    
    line(timesF, semiS, 'color', [1,0,0],'linewidth',3); hold off;     % Plot graph
    line(timesF, mainS, 'color', [0,0,0],'linewidth',3); hold on;      % Plot graph
    
    ylim([min(mainScores), max(mainScores)]+[-1,1]*0.5); xlim([0, max(timesF)]);
    t=1000*3600; set(gca, 'xtick', t:t:10e9, 'ytick', mainScores);
    set(gca, 'xticklabel', 1:24, 'yticklabel', titles);
    

end

function plotPercentageState(scoring, states, statesCode, stateColors)
    scr = zeros(1,length(states));
    for state_i_plot_Per = 1:length(statesCode)
        scr(state_i_plot_Per) = sum(scoring == statesCode(state_i_plot_Per))/length(scoring) * 100;
    end
    
    
    % Plot
    subplot(plotRows,plotCols,subplots{plotBars}); hold on
    for C = 1:length(scr)
        bar(C,scr(C),'FaceColor',stateColors(C,:));
    end
    hold off;  
end
        
function designSubplots()
    subplot(plotRows,plotCols,subplots{plotFFT1}); xlim([0 40]); xlabel('freqency'); ylabel('log power'); set(gca, 'xtick',0:10:40); title ('Frontal EEG','unit','normalized','Position',[0.5,0.85]);
    subplot(plotRows,plotCols,subplots{plotFFT2}); xlim([0 40]); xlabel('Freqency','unit','normalized','Position',[0.5,-0.016]); set(gca, 'xtick',0:10:40,'xticklabel',{'0','10','','30','40'}); ylabel('log power'); title ('Parietal EEG','unit','normalized','Position',[0.5,0.85]);
    subplot(plotRows,plotCols,subplots{plotScatter}); set(gca,'XColor',get(gcf,'Color'),'YColor',get(gcf,'Color'),'TickDir','out'); CHECK_YLIM = get(gca,'ylim');
    subplot(plotRows,plotCols,subplots{plotEEG}); xlabel('High / low EEG power ratio (log scale)'); ylabel('Occurence (%)');
    subplot(plotRows,plotCols,subplots{plotEMG}); set(gca,'view',[90 -90],'xlim',CHECK_YLIM); xlabel('EMG RMS (log scale)'); ylabel('Occurence (%)');
    subplot(plotRows,plotCols,subplots{plotBars}); set(gca, 'ytick',0:10:100); xlabel('States'); ylabel('Occurence (%)');
end

function fullScreen(f)
    set(f,'units','normalized','outerposition',[0 0 1 1]);
end

end