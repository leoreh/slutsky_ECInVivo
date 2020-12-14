function fepsp = fEPSP_analysis(varargin)

% gets a fepsp struct and analyzes the traces according to protocol
% (currently io or stp, in the future maybe more). assumes fEPSPfromDat or
% WCP has been called beforehand.
%
% INPUT
%   fepsp       struct. see fEPSPfromDat or fEPSPfromWCP
%   basepath    string. path to .dat file (not including dat file itself)
%   dt          numeric. deadtime for exluding stim artifact
%   force       logical. force reload {false}.
%   saveVar     logical. save variable {1}.
%   saveFig     logical. save graphics {1}.
%   graphics    numeric. if 0 will not plot grpahics. if greater than
%               nspkgrp will plot all grps, else will plot only
%               selected grp {1000}.
%   vis         char. figure visible {'on'} or not ('off')
%   savename    char. name to save the variable by. By default will be
%               folder name from 'basepath' + .fepsp {[]}.
%   MainTets    Numeric vec, Main Tetrodes, choose waveform edge points only for them.
%               All other Tetrodes will use an avergage of the times
%               choosen for them, in the matching intensity. If empty or 
%               any > fepsp.info.spkgrp, will treat all Tetrodes as important {[]}.
%
% CALLS
%   none
%
% OUTPUT
%   fepsp       struct with fields described below
%
% TO DO LIST
%   # Lior Da Marcas take over (done)
%   # Better warning msg to "Polynomial is not unique; degree >=number of data points"
%   # Add ablity to cheack diffrent time using LineChooseWin before saving?
%   # Add option to remove traces using LineChooseWin?
%   # Add option to re-analyse only some tetrodes?
%
% 16 oct 20 LH  UPDATES
% 02 Nov 20 LD  Change analysis from range on relevant window to
%               amplitude as 1st peak on waveArg, Change graphic to work
%               with that
% 05 Dec 20 LD  Change analysis from auto to user define Points via mini
%               GUI, add slope analysis, minor change graphics to comply
%               with this changes and to export maximazed view of graph.
%               Give option to caller to determine file name at saving
% 09 Dec 20 LD  Added MainTets, and now user choose waveform edge points
%               for each intensity per tetrode

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addOptional(p, 'fepsp', []);
addOptional(p, 'basepath', pwd);
addOptional(p, 'dt', 2, @isnumeric);
addOptional(p, 'force', false, @islogical);
addOptional(p, 'saveVar', true, @islogical);
addOptional(p, 'saveFig', true, @islogical);
addOptional(p, 'graphics', 1000);
addOptional(p, 'vis', 'on', @ischar);
addOptional(p, 'savename', [], @ischar);
addOptional(p, 'MainTets', [], @isnumeric);

parse(p, varargin{:})
fepsp = p.Results.fepsp;
basepath = p.Results.basepath;
dt = p.Results.dt;
force = p.Results.force;
saveVar = p.Results.saveVar;
saveFig = p.Results.saveFig;
graphics = p.Results.graphics;
vis = p.Results.vis;
savename = p.Results.savename;
MainTets = p.Results.MainTets;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get params from fepsp struct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~, basename] = fileparts(basepath);
if isempty(savename)
    fepspname = [basename '.fepsp.mat'];
else
    fepspname = savename;
end

% try to load file if not in input
if isempty(fepsp)
    if exist(fepspname, 'file')
        load(fepspname)
    end
end
if isfield(fepsp, 'slope_10_50') && ~force
    load(fepspname)
    return
end

fs = fepsp.info.fs;
spkgrp = fepsp.info.spkgrp;
nspkgrp = length(spkgrp);
nfiles = length(fepsp.intens);
protocol = fepsp.info.protocol;
% make sure tstamps column vector
fepsp.tstamps = fepsp.tstamps(:);
tstamps = fepsp.tstamps;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prepare for analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dt = round(dt / 1000 * fs);
switch protocol
    case 'io'
        % single pulse of 500 us after 30 ms. recording length 150 ms.
        % repeated once every 15 s. negative peak of response typically
        % 10 ms after stim.
        nstim = 1;
        [~, wvwin(1)] = min(abs(tstamps - 0));
        [~, wvwin(2)] = min(abs(tstamps - 30));
        wvwin(1) = wvwin(1) + dt;
    case 'stp'
        % 5 pulses of 500 us at 50 Hz. starts after 10 ms. recording length
        % 200 ms. repeated once every 30 s
        nstim = 5;
        
        % correct stim frequency
        if strcmp(fepsp.info.recSystem, 'oe') || strcmp(fepsp.info.recSystem, 'tdt')
                ts = fepsp.info.stimTs;
                ts = mean(ts(ts < 500)) / fs * 1000;
        elseif strcmp(fepsp.info.recSystem, 'wcp')
                ts = 20;               
        end     
        wvwin = round([10 : ts : 5 * ts; 30 : ts : 5 * ts + 10]' * fs / 1000);
        wvwin(:, 1) = wvwin(:, 1) + dt;
        wvwin(:, 2) = wvwin(:, 2) - dt;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% analyze
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% traceAvg      3d mat (tetrode x intensity x sample) of entire trace
fepsp.traceAvg  = nan(nspkgrp, nfiles, length(fepsp.tstamps));
% waves         2d cell (tetrode x intensity) where each cell contains the
%               waves (zoom in view of traces)
fepsp.waves     = cell(nspkgrp, nfiles);
% wavesAvg      3d mat (tetrode x intensity x sample) of waves (zoom in
%               view of trace), averages across traces. for io only
fepsp.wavesAvg  = nan(nspkgrp, nfiles, length(wvwin(1) : wvwin(2)));
% info.AnalysedTimePoints   3d mat (tetrode x intensity x stim), The time points user
%               choose for each stim. For each 2 of dim3  "columns", the first is the
%               start and the secound is the end of the area to analyse.
%               Time is from fepsp.tstamps.
fepsp.info.AnalysedTimePoints = nan(nspkgrp,nfiles,2*nstim);
fepsp.info.MainTets = MainTets;
% ampcell       2d array (tetrode x stim) where each cell contains the
%               amplitude/s for each trace
fepsp.ampcell   = cell(nspkgrp, nfiles);
% slopecell_10_50   2d array (tetrode x intensity) where each cell contains
%               the slope/s of 10% to 50% for each trace
fepsp.slopecell_10_50 = cell(nspkgrp, nfiles);
% slopecell_20_90   2d array (tetrode x intensity) where each cell contains
%               the slope/s of 20% to 90% for each trace
fepsp.slopecell_20_90 = cell(nspkgrp, nfiles);
% amp           2d (io) or 3d (stp) mat (tetrode x intensity x stim) of
%               amplitude averaged across traces
fepsp.amp   = nan(nspkgrp, nfiles, nstim);
% slope_10_50   2d (io) or 3d (stp) mat (tetrode x intensity x stim) of
%               slope 10% to 50% averaged across traces (STP not implanted!)
fepsp.slope_10_50 = nan(nspkgrp, nfiles, nstim);
% slope_20_90    2d (io) or 3d (stp) mat (tetrode x intensity x stim) of
%               slope 20% to 90% averaged across traces (STP not implanted!)
fepsp.slope_20_90 = nan(nspkgrp, nfiles, nstim);
% ampNorm       2d array of normalized amplitudes. for each trace the
%               responses are normalized to first response. these
%               normalized amplitudes. for
%               stp only
fepsp.ampNorm   = nan(nspkgrp, nfiles, nstim);
% facilitation  2d mat of average maximum normalized response. for stp only
fepsp.facilitation = nan(nspkgrp, nfiles);

fepsp = rmfield(fepsp, 'ampNorm');
The_10_50_20_90 = cell(nspkgrp,nfiles); %Slope ind for average waveforms (See later)
TimeFrameWindow = fepsp.tstamps(wvwin(1) : wvwin(2)); %Timeframe or responce 
StartStimTimeIND = nan(nspkgrp,nfiles);
EndStimTimeIND = nan(nspkgrp,nfiles);

% Choose waveform edge points for Main Tets, and take an avarage of it for all others
if isempty(MainTets) || any(MainTets > nspkgrp) % If MainTets is empty or any is too big, take all 
    MainTets = 1 : nspkgrp;
end
for j = MainTets
    yLimit = [min(fepsp.traces{j, end}, [], 'all'),...
        max(fepsp.traces{j, end}, [], 'all')];
    for i = nfiles : - 1 : 1
        % User choose timePoints for analyse by mini GUI, While all tet's traces are Presented
        %AllTraces = cell2mat(fepsp.traces(j,i));
        ChoosenTime = LineChooseWin(TimeFrameWindow,...
            fepsp.traces{j,i}(wvwin(1) : wvwin(2), :), nstim,...
            j, fepsp.intens(i), yLimit);
        
        % Calculate the closest Point that actually exist in data
        TimeFrameWindowMat = repmat(TimeFrameWindow, 1, nstim);
        [~,StartStimTimeIND(j, i)] = min(abs(TimeFrameWindowMat - ChoosenTime(1:2:end)), [], 1);
        [~,EndStimTimeIND(j, i)] = min(abs(TimeFrameWindowMat - ChoosenTime(2:2:end)), [], 1);
%         fepsp.info.AnalysedTimePoints(j,1:2:end) = TimeFrameWindow(StartStimTimeIND(j,i));
%         fepsp.info.AnalysedTimePoints(j,2:2:end) = TimeFrameWindow(EndStimTimeIND(j,i));
    end
end
DoneTetVec = ismember(1:nspkgrp,MainTets);
StartStimTimeIND(~DoneTetVec,:) = repmat(round(mean(StartStimTimeIND(DoneTetVec,:),1)),sum(~DoneTetVec),1);
EndStimTimeIND(~DoneTetVec,:) = repmat(round(mean(EndStimTimeIND(DoneTetVec,:),1)),sum(~DoneTetVec),1);
fepsp.info.AnalysedTimePoints(:,:,1:2:end) = TimeFrameWindow(StartStimTimeIND);
fepsp.info.AnalysedTimePoints(:,:,2:2:end) = TimeFrameWindow(EndStimTimeIND);

% Calculate by protocol
for j = 1 : nspkgrp
    for i = 1 : nfiles
        fepsp.traceAvg(j, i, :) = mean(fepsp.traces{j, i}, 2);
        switch protocol
            case 'io'
                fepsp.waves{j, i} = fepsp.traces{j, i}(wvwin(1) : wvwin(2), :);
                % Amp is simply absolute Start point - End point
                fepsp.ampcell{j, i} = abs(fepsp.waves{j, i}(StartStimTimeIND(j,i),:) - fepsp.waves{j, i}(EndStimTimeIND(j,i),:));
                
                % Calculate  Each trace slopes
                for kk = length(fepsp.ampcell{j, i}):-1:1
                    % Find for each trace the closest points that define about 10%-50% & 20%-90% of the calculated amp.
                    % Can be vectorize relatively simply using 3d array if run time too slow (Unexpected).
                    AllPointsAmp = abs(fepsp.waves{j, i}(StartStimTimeIND(j,i):EndStimTimeIND(j,i),kk) - fepsp.waves{j, i}(StartStimTimeIND(j,i),kk));
                    MinusMat = fepsp.ampcell{j, i}(kk)*[0.1 0.5 0.2 0.9];
                    [~,Trace_10_50_20_90] = min(abs(AllPointsAmp-MinusMat),[],1);
                    % If there are both a significant biphasic between StartStim & EndStim, the lower precentage might be found after
                    % the higher one. Will break colon operator later. in this case, just flip them.
                    if Trace_10_50_20_90(1) > Trace_10_50_20_90(2)
                        Trace_10_50_20_90(1:2) = Trace_10_50_20_90(2:-1:1);
                    end
                    if Trace_10_50_20_90(3) > Trace_10_50_20_90(4)
                        Trace_10_50_20_90(3:4) = Trace_10_50_20_90(4:-1:3);
                    end
                    Trace_10_50_20_90 = Trace_10_50_20_90+StartStimTimeIND(j,i)-1; %Match the indicies from small StartStim:EndStim window back to TimeFrameWindow.
                    
                    % Calculate Slope by fitting a line. Cannot be vectorize simply.
                    FitParams = polyfit(TimeFrameWindow(Trace_10_50_20_90(1) :...
                        Trace_10_50_20_90(2)),...
                        fepsp.waves{j, i}(Trace_10_50_20_90(1) :...
                        Trace_10_50_20_90(2), kk), 1);
                    fepsp.slopecell_10_50{j, i}(kk) = FitParams(1);
                    FitParams = polyfit(TimeFrameWindow(Trace_10_50_20_90(3) :...
                        Trace_10_50_20_90(4)),...
                        fepsp.waves{j, i}(Trace_10_50_20_90(3) :...
                        Trace_10_50_20_90(4), kk), 1);
                    fepsp.slopecell_20_90{j, i}(kk) = FitParams(1);
                end
                
                %Redo on mean wave. Save closest points that define about 10%-50% & 20%-90% to avoid recalcuating for graphics
                fepsp.wavesAvg(j, i, :) = mean(fepsp.waves{j, i}, 2);
                fepsp.amp(j, i) = abs(fepsp.wavesAvg(j, i, StartStimTimeIND(j,i)) - fepsp.wavesAvg(j, i,EndStimTimeIND(j,i)));
                
                AllPointsAmp = squeeze(abs(fepsp.wavesAvg(j, i, StartStimTimeIND(j,i):EndStimTimeIND(j,i)) - fepsp.wavesAvg(j, i, StartStimTimeIND(j,i))));
                MinusMat = fepsp.amp(j, i)*[0.1 0.5 0.2 0.9];
                [~,The_10_50_20_90{j,i}] = min(abs(AllPointsAmp-MinusMat),[],1);
                if The_10_50_20_90{j,i}(1) > The_10_50_20_90{j,i}(2)
                        The_10_50_20_90{j,i}(1:2) = The_10_50_20_90{j,i}(2:-1:1);
                end
                if The_10_50_20_90{j,i}(3) > The_10_50_20_90{j,i}(4)
                        The_10_50_20_90{j,i}(3:4) = The_10_50_20_90{j,i}(4:-1:3);
                end
                The_10_50_20_90{j,i} = The_10_50_20_90{j,i}+StartStimTimeIND(j,i)-1;
                FitParams = polyfit(TimeFrameWindow(The_10_50_20_90{j,i}(1) :...
                    The_10_50_20_90{j,i}(2)),...
                    squeeze(fepsp.wavesAvg(j, i, The_10_50_20_90{j, i}(1) :...
                    The_10_50_20_90{j,i}(2))), 1);
                fepsp.slope_10_50(j, i) = FitParams(1);
                FitParams = polyfit(TimeFrameWindow(The_10_50_20_90{j,i}(3) :...
                    The_10_50_20_90{j,i}(4)),...
                    squeeze(fepsp.wavesAvg(j, i, The_10_50_20_90{j, i}(3) :...
                    The_10_50_20_90{j,i}(4))),1);
                fepsp.slope_20_90(j, i) = FitParams(1);
            case 'stp'
                % note; after reviewing the data it seems that specifically
                % for stp maximum absolute value may be better than range
                for ii = 1 : nstim
                    fepsp.ampcell{j, i}(ii, :) =...
                        range(fepsp.traces{j, i}(wvwin(ii, 1) :  wvwin(ii, 2), :));
                end
                fepsp.ampNorm{j, i} = fepsp.ampcell{j, i} ./ fepsp.ampcell{j, i}(1, :);
                fepsp.facilitation(j, i) = mean(max(fepsp.ampNorm{j, i}));
        end
    end
end

% save updated struct
if saveVar
    save(fepspname, 'fepsp');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if graphics
    if graphics > nspkgrp
        grp = 1 : nspkgrp;
    else
        grp = graphics;
    end
    for i = grp
        switch protocol
            case 'io'
                fh = figure('Visible', vis,'WindowState','maximized'); %Max window for better looking export
                suptitle(sprintf('T#%d', i)) 
                subplot(3, 1, 1) %Avarge Waveforms Marked
                plot(TimeFrameWindow, squeeze(fepsp.wavesAvg(i, :, :))','LineWidth',1)
                hold on
                % Mark on each waveform Start & End of analysed area
                Yind = sub2ind(size(fepsp.wavesAvg),ones(1,size(fepsp.wavesAvg,2))*i,1:size(fepsp.wavesAvg,2),StartStimTimeIND(i,:));
                plot(fepsp.info.AnalysedTimePoints(i,:,1:2:end),fepsp.wavesAvg(Yind),'*')
                Yind = sub2ind(size(fepsp.wavesAvg),ones(1,size(fepsp.wavesAvg,2))*i,1:size(fepsp.wavesAvg,2),EndStimTimeIND(i,:));
                plot(fepsp.info.AnalysedTimePoints(i,:,2:2:end),fepsp.wavesAvg(Yind),'*')
                % Marker area in which slope was analysed
                for jj = 1:length(The_10_50_20_90(i,:))
                    P_10_50 = plot(TimeFrameWindow(The_10_50_20_90{i,jj}(1):The_10_50_20_90{i,jj}(2)),...
                        squeeze(fepsp.wavesAvg(i,jj,The_10_50_20_90{i,jj}(1):The_10_50_20_90{i,jj}(2))),'b','LineWidth',3);
                    P_10_50.Color(4) = 0.5;
                    P_10_50.Annotation.LegendInformation.IconDisplayStyle = 'off';
                    P_20_90 = plot(TimeFrameWindow(The_10_50_20_90{i,jj}(3):The_10_50_20_90{i,jj}(4)),...
                        squeeze(fepsp.wavesAvg(i,jj,The_10_50_20_90{i,jj}(3):The_10_50_20_90{i,jj}(4))),'y','LineWidth',3);
                    P_20_90.Color(4) = 0.5;
                    P_20_90.Annotation.LegendInformation.IconDisplayStyle = 'off';
                end
                axis tight
                ylim(ylim()*1.05)
%                 yLimit = [min([fepsp.wavesAvg(:)]) max([fepsp.wavesAvg(:)])];
%                 ylim(yLimit)
                xlabel('Time [ms]')
                ylabel('Voltage [mV]')
                legend([split(num2str(sort(fepsp.intens)));{'Amp measure point Stim Start'};{'Amp measure point Stim End'}],'Location','best','NumColumns',2)
                box off
                P_10_50.Annotation.LegendInformation.IconDisplayStyle = 'on';
                P_10_50.DisplayName = 'Slope Area 10%-50%';
                P_20_90.Annotation.LegendInformation.IconDisplayStyle = 'on';
                P_20_90.DisplayName = 'Slope Area 20%-90%';
                
                subplot(3, 1, 2) %Amplitude across ints (on avarge waveforms)
                %ampmat = cell2nanmat(fepsp.ampcell(i, :));
                %boxplot(ampmat, 'PlotStyle', 'traditional')
                bar(fepsp.amp(i,:))
                ylim([min(horzcat(fepsp.ampcell{:})) max(horzcat(fepsp.ampcell{:}))])
                xticklabels(split(num2str(sort(fepsp.intens))))
                xlabel('Intensity [uA]')
                ylabel('Amplidute [mV]')
                box off
                
                subplot(3, 1, 3) %Slope across ints (on avarge waveforms)
                bar([fepsp.slope_10_50(i,:);fepsp.slope_20_90(i,:)])
                xticklabels({'From 10% to 50%' 'From 20% to 90%'})
                xlabel('Measure area')
                ylabel('Slope [mV/mS]')
                legend([split(num2str(sort(fepsp.intens)))])
                box off
                
            case 'stp'
                fh = figure('Visible', vis);
                suptitle(sprintf('%s - T#%d', basename, i))
                subplot(2, 1, 1)
                plot(fepsp.tstamps, squeeze(fepsp.traceAvg(i, :, :))')
                axis tight
                yLimit = [min(min(horzcat(fepsp.traces{i, :})))...
                    max(max(horzcat(fepsp.traces{i, :})))];
                ylim(yLimit)
                hold on
                plot(repmat([0 : ts : ts * 4]', 1, 2), yLimit, '--k')
                xlabel('Time [ms]')
                ylabel('Voltage [mV]')
                legend(split(num2str(sort(fepsp.intens))))
                box off
                
                subplot(2, 1, 2)
                for ii = 1 : length(fepsp.intens)
                    x(ii, :) = mean(fepsp.ampNorm{i, ii}, 2);
                end
                plot([1 : nstim], x)
                xticks([1 : nstim])
                xlabel('Stim No.')
                ylabel('Norm. Amplitude')
                yLimit = ylim;
                ylim([0 yLimit(2)])
        end
        if saveFig
            figpath = fullfile(basepath, 'graphics');
            mkdir(figpath)
            figname = [figpath '\' basename '_fepsp_t' num2str(i)];
            export_fig(figname, '-tif', '-r300', '-transparent')
        end
    end
end

end

%% Helper LineChooseWin function
function [LineXPos] = LineChooseWin(Data2PlotX,Data2PlotY,nstim,TetNum,Intens,yLimit)
%[LineXPos] = LineChooseWin(Data2PlotX,Data2PlotY,nstim,TetNum)
%Helper function to fEPSP_analysis, Plot data and choose stimuli time
%   Will plot the inserted data and let you choose start & and time for
%   each stimuli by moving lines ob the plot.
%   Inputs:
%       Data2PlotX - numeric, data for x axis (usally double time). See "plot" for limits.
%       Data2PlotY - numeric, data for y axis. See "plot" for limits.
%       nstim      - numeric scalar, number of stimuli to choose there time.
%       TetNum     - numeric scalar, the number of tetrode corrently working on,
%                    just for title.
%   Outputs:
%       LineXPos   - numeric vector, lenght 2*nstim. the x position of each
%                    line in refernce to Data2PlotX. The Vector is ordered
%                    StartStim1 EndStim1 StartStim2 EndStim2 and so on.

% Main Graph
ChooseWin = figure('WindowState','maximized');
plot(Data2PlotX,Data2PlotY);
xlabel('Time [ms]')
ylabel('Amp [mV]')
title({sprintf('T%d - %d uA', TetNum, Intens)...
    'Move green / red lines to start / end of each response'...
 'Press "Enter\return" when done'})
axis tight

% Creating the Lines
XLims = xlim();
ylim(yLimit);
ylim(ylim()*1.1);
LineStartPos = linspace(XLims(1),XLims(2),2*nstim);
Colors = {'g','r'};
Pos = {'S','E'};
for ii = length(LineStartPos):-1:1
    ColorIND = mod(ii,2);
    ColorIND(ColorIND == 0) = 2;
    d(ii) = drawline('Position',[repmat(LineStartPos(ii),2,1) ylim()'],'InteractionsAllowed','translate','Color',Colors{ColorIND},...
        'Label',sprintf('%d%s t=%.1f',ceil(ii/2),Pos{ColorIND},LineStartPos(ii)));
end
Lis = addlistener(d,'MovingROI',@WriteLoc); %Listener for showing the time in lines lable

%Waiting to user to place lines, and making comfy close options, all simply
%continue the function, to prevent closing before saving lines pos 
uicontrol(ChooseWin,'Style','pushbutton','String','Finished','Callback',@(~,~) uiresume)
ChooseWin.CloseRequestFcn = @(~,~) uiresume; %Normal close will not return pos, so just use the uiresume when closing
ChooseWin.KeyReleaseFcn = @ResumeIfEnter;
uiwait(ChooseWin);

%Saving lines Xpos then close the figure
LinePos = [d.Position];
LineXPos = LinePos(1,1:2:end);
delete(Lis)
delete(ChooseWin)

    function WriteLoc(obj,evt)
        %Callback function to the listner, just change the tie written in
        %line lable when moving it
        CurrXPos = evt.CurrentPosition(1);
        obj.Label = [obj.Label(1:(find(obj.Label == '='))), sprintf('%.1f',CurrXPos)];
    end
    function ResumeIfEnter(~,evt)
        %Figure KeyReleaseFcn, just make pressing "enter"\"return" resume
        %the function (therefore saving lines XPos and closing figure)
        if evt.Character == 13
            uiresume
        end
    end
end

% EOF