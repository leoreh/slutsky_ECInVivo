function ied = analyze(ied,varargin)
% analyze the ied result for rate, and display figure with info about the ied
%
%   INPUT (in 1st position):
%       IED.data object, after detection.
%   INPUT (optional, name value):
%       binsize     scalar {60*ied.fs} in [samples]. for rate calculation.
%       smf         smooth factor for rate [bins] {7}.
%       marg        scalar {0.05} in [s]. time margin for clipping discharges.
%       saveVar     logical {true}. save variable.
%       basepath    recording session path {pwd}, folder to save in.
%       basename    string. if empty extracted from basepath. File name to save in.
%       saveFig     logical {true}. save figure
% **IMPORTANT:** When a parameter exists in both `IED.data` and is also
%                provided by the user as a name-value pair, the user's
%                parameter value will take precedence, overwriting the
%                corresponding value in `IED.data`.
%
% OUTPUT
%   IED.data object
%
% DEPENDECIES: (out of package)
%   slutsky_ECInVivo\utilities\times2rate.m
%   slutsky_ECInVivo\spikes\correlation\CCG.m
%   export_fig (not in slutsky_ECInVivo, only needed for support of 2020a<)
%
% Based on getIIS by LH (see +IED/legacy folder)
% By: LdM
% Published: 230827
%
%   see also IED.data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isa(ied,'IED.data')
    error("first input must be IED.data obj")
end

p = inputParser;
addParameter(p, 'binsize', 60*ied.fs, @isnumeric)
addParameter(p, 'marg', 0.05, @isnumeric)
addParameter(p, 'smf', 7, @isnumeric)
addParameter(p, 'saveVar', true, @islogical);
addParameter(p, 'basepath', pwd, @isstr);
addParameter(p, 'basename', [], @isstr);
addParameter(p, 'saveFig', true, @islogical);
parse(p, varargin{:})


if  ~ismember("binsize",p.UsingDefaults) || isempty(ied.binsize)
    % user ask for binsize, or no specific binsize exist in ied -
    % overwrite existing
    ied.binsize = p.Results.binsize;
end
if ~ismember("smf",p.UsingDefaults) || isempty(ied.smf)
    % user ask for smf, or no specific smf exist in ied -
    % overwrite existing
    ied.smf = p.Results.smf;
end
if ~ismember("marg",p.UsingDefaults) || isempty(ied.marg)
    % user ask for marg, or no specific marg exist in ied -
    % overwrite existing
    ied.marg = p.Results.marg;
end

saveVar = p.Results.saveVar;
basepath = p.Results.basepath;
basename = p.Results.basename;
saveFig = p.Results.saveFig;

if isempty(basename)
    % write where to save
    [~,basename] = fileparts(basepath);
end
if any(~ismember(["basepath","basename"],p.UsingDefaults)) && (saveVar || saveFig)
    % if user gave any info & want to save, overwrite existing
    ied.file_loc = fullfile(basepath,join([basename "ied.mat"],"."));
end

% warn if ied wasn't curated
if ~ismember(ied.status,["curated", "analyzed", "changed_post_analysis"])
    warning('ied wasn''t manually curated before analysis! Run IED.curate to perform.')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% remove discharges with bad margins
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

margs = floor(ied.marg * ied.fs);           % margs [samples]; ied.marg [ms]
accepted = ied.accepted;
out_margs = find( (ied.pos+margs > length(ied.sig)) | (ied.pos-margs < 1) );
if ~isempty(out_margs)
    accepted(out_margs) = false;
    fprintf("\n**** Making discharges {%s} not-accepted for analysis, as their margins is partly out of signal ****\n",...
        num2str(out_margs,"%d, "))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calc rate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% collect discharges that passed curation
true_pos = ied.pos(accepted);

[ied.rate, ied.edges, ied.cents] = times2rate(true_pos, 'winCalc', [1, length(ied.sig)],...
    'binsize', ied.binsize, 'c2r', false);
% this is super dangerous because if binsize < 1 min than the rate will
% effectively be greater than the number of counts
% iis.rate = iis.rate * fs * 60;      % convert counts in bins to 1 / min
ied.rate = movmean(ied.rate, ied.smf);
ied.status = "analyzed";

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save result if needed (overwrite)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if saveVar
    save(ied.file_loc,"ied")
    fprintf("\n****** Save in %s ******\n",ied.file_loc)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prepare for graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% collect clips for displaying later
sig_tstamps =  (1 : length(ied.sig))' / ied.fs;
[clipped_discharges, clipped_tstamps] = extract_discharges(ied);
peak_val = ied.sig(true_pos);

%%%% prepare power spectrum via wavelet
% average wavelet coefficients for each discharge. note this produces very
% similar results to the coefficients obtained from the average discharge
% waveform. individual coefficients may still be carried out if the width
% of each discharge is to be calculated

fb = cwtfilterbank('SignalLength', size(clipped_discharges, 2), 'VoicesPerOctave', 32,...
    'SamplingFrequency', ied.fs, 'FrequencyLimits', [1 ied.fs / 4]);
mwv = mean(clipped_discharges, 1)'; % mean waveform
[cfs, f, coi] = cwt(mwv, 'FilterBank', fb);

%%%% prepare autocorralations (ACG)
[ccg, tccg] = CCG({true_pos / ied.fs}, [], 'duration', 10,...
    'binSize', 0.1);
[ccg2, tccg2] = CCG({true_pos / ied.fs}, [], 'duration', 30,...
    'binSize', 0.3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fh = figure('Color','w','UserData',ied);
fh.UserData = struct('ied',ied,'confirm_decline',true);
set(fh, 'units','normalized','outerposition',[0 0 1 1]);
sgtitle(basename,'Interpreter','none')

%%%% Prepare uicontextmenus
cm_enlarge = uicontextmenu(fh);
cm_decline_discharge = uicontextmenu(fh);
uimenu(cm_enlarge,"Text","Enlarge","MenuSelectedFcn",@undock_axe)
uimenu(cm_decline_discharge,"Text","Decline Event","MenuSelectedFcn",@decline_discharge)

% raw and ied rate
ax = subplot(3, 4, 1 : 2);
plot(sig_tstamps / 60, ied.sig ,'k')
axis tight
hold on
markthr(ax,ied.thrDir,ied.thr(2))
ylabel('Voltage [mV]')
yyaxis right
plot(ied.cents / ied.fs / 60, ied.rate, 'b', 'LineWidth', 3)
ylabel('Rate [IEDs / bin]')
set(gca, 'TickLength', [0 0])
box off
title('Raw signal and IIS rate')
ax.UIContextMenu = cm_enlarge;

% detection
ax = subplot(3, 4, 5 : 6);
plot(sig_tstamps / 60, zscore(ied.sig), 'k')
hold on
axis tight
markthr(ax,ied.thrDir,ied.thr(1))
plot([true_pos true_pos] / ied.fs / 60, [-10 -1], '--g', 'LineWidth', 2)
xlabel('Time [m]')
ylabel('Z-score')
set(gca, 'TickLength', [0 0])
box off
title('IED detection')
ax.UIContextMenu = cm_enlarge;

% zoom in
ax = subplot(3, 4, 9 : 10);
midsig = round(length(ied.sig) / 2);
idx = round(midsig - 2 * ied.fs * 60 : midsig + 2 * ied.fs * 60);
idx2 = true_pos > idx(1) & true_pos < idx(end);
plot(sig_tstamps(idx) / 60, ied.sig(idx), 'k')
axis tight
hold on
scatter(true_pos(idx2) / ied.fs / 60,...
    peak_val(idx2), '*');
ylabel('Voltage [mV]')
xlabel('Time [m]')
xticks(round([midsig / ied.fs / 60 - 2, midsig / ied.fs / 60 + 2]))
set(gca, 'TickLength', [0 0])
box off
title('Mid Signal Zoom In')
ax.UIContextMenu = cm_enlarge;

% IED waveforms
ax = subplot(3, 4, 3);
p = plot(clipped_tstamps * 1000, clipped_discharges');
[~,pos_idx] = ismember(true_pos,ied.pos);
for iLine = 1:numel(p)
    p(iLine).UIContextMenu = cm_decline_discharge;
    p(iLine).Tag = string(pos_idx(iLine));
    p(iLine).DataTipTemplate.DataTipRows(end+1) = ...
        dataTipTextRow('IED:',repelem(string(pos_idx(iLine)),numel(p(iLine).XData)));
end
pos_idx = num2cell(string(pos_idx));
[p.Tag] = deal(pos_idx{:});
ylabel('Voltage [mV]')
xlabel('Time [ms]')
axis tight
xticks([-ied.marg, 0, ied.marg] * 1000);
set(gca, 'TickLength', [0 0])
box off
title('IED waveform')
ax.UIContextMenu = cm_enlarge;

% mean + std waveform
ax = axes('Position',[.542 .71 .09 .07]);
box on
IED.utils.stdshade(ax,double(clipped_discharges), 0.5, 'k', clipped_tstamps);
axis tight
set(gca, 'TickLength', [0 0], 'YTickLabel', [], 'XTickLabel', [],...
    'XColor', 'none', 'YColor', 'none', 'Color', 'none')
title(sprintf('n = %d', size(clipped_discharges, 1)));
box off
ax.UIContextMenu = cm_enlarge;

% ied cwt
ax = subplot(3, 4, 4);
imagesc(clipped_tstamps, f, abs(cfs))%,'HitTest','off')
axis xy
set(gca, 'YScale', 'log')
% caxis([0.0001 1])
origSize = get(gca, 'Position');
colorbar
set(gca, 'Position', origSize);
hold on
plot(clipped_tstamps, coi, 'w', 'LineWidth', 2)
ylim([min(f) max(f)]);
xlim([clipped_tstamps(1) clipped_tstamps(end)]);
yticks([ceil(min(f)), 100, max(f)])
xticks([min(clipped_tstamps), 0, max(clipped_tstamps)])
ylabel('Frequency [Hz]')
xlabel('Time [ms]')
set(gca, 'TickLength', [0 0])
box off
title('Average Scalogram')
ax.UIContextMenu = cm_enlarge;

% amplitude histogram
ax = subplot(3, 4, 8);
h = histogram(log10(abs(double(peak_val))), 30, 'Normalization', 'Probability');
h.EdgeColor = 'none';
h.FaceColor = 'k';
h.FaceAlpha = 1;
xlabel('Peak voltage [log(uV)]')
ylabel('Probability [%]')
set(gca, 'TickLength', [0 0])
box off
title('Amplitude Distribution')
ax.UIContextMenu = cm_enlarge;

% max frequency and amplitude vs. time
ax = subplot(3, 4, 7);
scatter(true_pos / ied.fs / 60, double(peak_val), 2, 'b', 'filled');
axis tight
x = xlim;
l2 = lsline;
set(l2, 'color', 'b')
l2.LineWidth = 3;
xlim(x);
ylabel('Amplitude [mV]')
axis tight
xlabel('Time [m]')
set(gca, 'TickLength', [0 0])
box off
title('Amplitude. vs. Time')
ax.UIContextMenu = cm_enlarge;

% ACH 5s
ax = subplot(3, 4, 11);
IED.utils.plotCCG('ccg', ccg, 't', tccg / 1000, 'basepath', basepath,...
    'saveFig', false, 'c', {'k'});
xlabel('Time [s]')
title('Autocorrelogram')
ax.UIContextMenu = cm_enlarge;

% ACH 10s
ax = subplot(3, 4, 12);
IED.utils.plotCCG('ccg', ccg2, 't', tccg2 / 1000, 'basepath', basepath,...
    'saveFig', false, 'c', {'k'});
xlabel('Time [s]')
title('Autocorrelogram')
ax.UIContextMenu = cm_enlarge;

if saveFig
    figname = replace(ied.file_loc,'.ied.mat','_IED');
    if exist("exportgraphics","file")
        figname = join([figname '.tif'],'');
        exportgraphics(fh,figname,"BackgroundColor","none")
    else
        export_fig(figname, '-tif', '-transparent')
    end

    %         savePdf(figname, basepath, fh)
end

end

function markthr(ax,thrDir,thr)
% simply add threshold marking to requested ax.
% INPUTS:
%   ax     - axis 2 mark on.
%   thrDir - string scalar, dircation of threshold, "positive",negative or "both".
%   thr    - threshold value.



% add threshold markers
if ismember(thrDir,["positive","both"])
    yline(ax, thr, '--r');
end
if ismember(thrDir,["negative","both"])
    yline(ax, -thr, '--r');
end
end

function undock_axe(src,~)
% create a new figure that is an enlargment of the existing axis.

%%%% collect which obj to work on
fig = src.Parent.Parent; % src = uimenue > uicontectmenu > figure
ax2enlarge = fig.CurrentObject;

%%%% copy obj to new figure
new_fig = figure();
try
    copied_ax = copyobj(ax2enlarge,new_fig,'legacy');
catch err
    if strcmp(err.message,'Object Copy of Axes with multiple coordinate systems is not supported.')
        % copy will fail when using yyaxis. Instead, create axe, copy each obj seperetly
        copied_ax = axes(new_fig);

        % orginaze left axis
        yyaxis(ax2enlarge,"left")
        %             yyaxis(copied_ax,"left") % no need - it is already on left
        copyobj(ax2enlarge.Children,copied_ax);
        %             copyobj(ax2enlarge.YLabel,copied_ax); % sometimes create bad label placment. just use the same text
        ylabel(copied_ax,ax2enlarge.YLabel.String)
        copied_ax.YLim = ax2enlarge.YLim;

        % orginaze tight axis
        yyaxis(ax2enlarge,"right")
        yyaxis(copied_ax,"right")
        copied_ax.YAxis(2).Color = ax2enlarge.YAxis(2).Color;
        copyobj(ax2enlarge.Children,copied_ax);
        %             copyobj(ax2enlarge.YLabel,copied_ax); % sometimes create bad label placment. just use the same text
        ylabel(copied_ax,ax2enlarge.YLabel.String)
        copied_ax.YLim = ax2enlarge.YLim;
        xlim(copied_ax,ax2enlarge.XLim)
    end
end

%%%% final touches

% make axis full size & look nice
copied_ax.Position = [0.1300 0.1100 0.7750 0.8150];
copied_ax.Color = 'w';

% restore context menu. assume that all ax children objects has 1 (or non has), and it is the same 1
cm2add = get(ax2enlarge.Children,'UIContextMenu');
if iscell(cm2add)
    cm2add = [cm2add{:}];
end
cm2add = unique(cm2add);
if ~isempty(cm2add)
    new_cm = copyobj(cm2add,new_fig,"legacy");
    [copied_ax.Children.UIContextMenu] = deal(new_cm);
end
% copy userdata for newfig context menu managers
new_fig.UserData = fig.UserData;
new_fig.UserData.OrigFig = fig;

% copy datatips
for iChild = numel(copied_ax.Children):-1:1
    if isprop(copied_ax.Children(iChild),'DataTipTemplate')
        copied_ax.Children(iChild).DataTipTemplate.DataTipRows = ax2enlarge.Children(iChild).DataTipTemplate.DataTipRows;
    end
end

end

function decline_discharge(src,~)
% mark a specific discharge as decline in original IED.data obj, and delete
% it from display.

% collect which obj to work on 
fig = src.Parent.Parent; % src = uimenue > uicontectmenu > figure
line2decline = fig.CurrentObject;
idx2change = str2double(line2decline.Tag);

% validate that user is sure about declination
if fig.UserData.confirm_decline
    answer = questDlg(...
        ["Are you sure you want to decline this spike?",...
        "This will affect the original IED.data!",...
        "Run IED.analyze again to calculate the change in rate & all"],...
        "Decline IED","Yes","Yes, do not ask again","No","No");
else
    % user ask not to ask again - always assume yes.
    answer = 'Yes';
end

if answer == "Yes, do not ask again"
    % Mark user request, always just decline from now
    fig.UserData.confirm_decline = false;
    if isfield(fig.UserData,"OrigFig")
        fig.UserData.OrigFig.UserData.confirm_decline = false;
    end
end
if contains(answer,"Yes")
    % decline discharge
    
    % change ied info - it isn't analyzed anymore
    ied = fig.UserData.ied;
    ied.accepted(idx2change) = false;
    [ied.rate, ied.edges, ied.cents] = deal([]);
    ied.status = "changed_post_analysis";
    
    % delete lines from display
    delete(line2decline);
    if isfield(fig.UserData,"OrigFig")
        line2decline = findobj(fig.UserData.OrigFig,"Tag",num2str(idx2change));
        delete(line2decline)
    end
end
end

% EOF