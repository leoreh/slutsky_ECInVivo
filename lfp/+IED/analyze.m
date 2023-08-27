function ied = analyze(ied,varargin)
% analyse the ied result for rate, and display figure with info about the ied
%
%   INPUT (in 1st position):
%       IED.data object, after detection.
%   INPUT (optional, name value):
%       smf         smooth factor for rate [bins] {7}.
%       saveVar     logical {true}. save variable
%       basepath    recording session path {pwd}
%       basename    string. if empty extracted from basepath
%       saveFig     logical {true}. save figure
%
% OUTPUT
%   IED.data object
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
addParameter(p, 'smf', 7, @isnumeric)
addParameter(p, 'saveVar', true, @islogical);
addParameter(p, 'basepath', pwd, @isstr);
addParameter(p, 'basename', [], @isstr);
addParameter(p, 'saveFig', true, @islogical);
parse(p, varargin{:})

ied.smf = p.Results.smf;
saveVar = p.Results.saveVar;
basepath = p.Results.basepath;
basename = p.Results.basename;
saveFig = p.Results.saveFig;

if (saveVar || saveFig)
    % write where to save
    if isempty(ied.file_loc)
        % create from user input
        if isempty(basename)
            [~,basename] = fileparts(basepath);
        end
        ied.file_loc = fullfile(basepath,join([basename "ied.mat"],"."));

    elseif ~isempty(ied.file_loc) && any(~ismember(["basename","basepath"],p.UsingDefaults))
        % inform user that its input did not matter
        fprintf("\n****Using already existing file path, as exist in ied (1st) input.\n" + ...
            "To save in a new path, change path before calling analyze function.")
    end


end

% warn if ied wasn't curated
if ~ismember(ied.status,["curated", "analysed"])
    warning('ied wasn''t manualy curated before analysis! Run IED.curate to perform.')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calc rate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% collect spikes that passed curation
true_pos = ied.pos(ied.accepted);

[ied.rate, ied.edges, ied.cents] = times2rate(true_pos, 'winCalc', [1, length(ied.sig)],...
    'binsize', ied.binsize, 'c2r', false);
% this is super dangerous because if binsize < 1 min than the rate will
% effectively be greater than the number of counts
% iis.rate = iis.rate * fs * 60;      % convert counts in bins to 1 / min
ied.rate = movmean(ied.rate, ied.smf);
ied.status = "analysed";
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save result if needed (overwrite)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if saveVar
    save(ied.file_loc,"ied")
    fprintf("\n****** Save in %s ******\n",ied.file_loc)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prapare for graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% collect clips for displaying later
sig_tstamps =  (1 : length(ied.sig))' / ied.fs;
[clipped_discharges, clipped_tstamps] = extract_discharges(ied);
peak_val = ied.sig(true_pos);

%%%% Prapare power spectrum via wavelet
% average wavelet coefficients for each spike. note this produces very
% similar results to the coefficients obtained from the average spike
% waveform. individual coefficients may still be carried out if the width
% of each spike is to be calculated

fb = cwtfilterbank('SignalLength', size(clipped_discharges, 2), 'VoicesPerOctave', 32,...
    'SamplingFrequency', ied.fs, 'FrequencyLimits', [1 ied.fs / 4]);
mwv = mean(clipped_discharges, 1)'; % mean waveform
[cfs, f, coi] = cwt(mwv, 'FilterBank', fb);

%%%% Prapare autocorralations (ACG)
[ccg, tccg] = CCG({true_pos / ied.fs}, [], 'duration', 10,...
    'binSize', 0.1);
[ccg2, tccg2] = CCG({true_pos / ied.fs}, [], 'duration', 30,...
    'binSize', 0.3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fh = figure('Color','w');
set(fh, 'units','normalized','outerposition',[0 0 1 1]);
sgtitle(basename)

% raw and iis
ax = subplot(3, 4, 1 : 2);
plot(sig_tstamps / 60, ied.sig ,'k')
axis tight
hold on
markthr(ied, ax)
ylabel('Voltage [mV]')
yyaxis right
plot(ied.cents / ied.fs / 60, ied.rate, 'b', 'LineWidth', 3)
ylabel('Rate [spikes / bin]')
set(gca, 'TickLength', [0 0])
box off
title('Raw signal and IIS rate')

% detection
ax = subplot(3, 4, 5 : 6);
plot(sig_tstamps / 60, zscore(ied.sig), 'k')
hold on
axis tight
markthr(ied, ax)
plot([true_pos true_pos] / ied.fs / 60, [-10 -1], '--g', 'LineWidth', 2)
xlabel('Time [m]')
ylabel('Z-score')
set(gca, 'TickLength', [0 0])
box off
title('IIS detection')

% zoom in
subplot(3, 4, 9 : 10)
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
title('IIS')

% iis waveforms
subplot(3, 4, 3)
plot(clipped_tstamps * 1000, clipped_discharges')
ylabel('Voltage [mV]')
xlabel('Time [ms]')
axis tight
xticks([-ied.marg, 0, ied.marg] * 1000);
set(gca, 'TickLength', [0 0])
box off
title('Spike waveform')

% mean + std waveform
axes('Position',[.542 .71 .09 .07])
box on
stdshade(clipped_discharges, 0.5, 'k', clipped_tstamps)
axis tight
set(gca, 'TickLength', [0 0], 'YTickLabel', [], 'XTickLabel', [],...
    'XColor', 'none', 'YColor', 'none', 'Color', 'none')
title(sprintf('n = %d', size(clipped_discharges, 1)));
box off

% iis cwt
subplot(3, 4, 4)
imagesc(clipped_tstamps, f, abs(cfs))
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

% amplitude histogram
subplot(3, 4, 8)
h = histogram(log10(abs(peak_val)), 30, 'Normalization', 'Probability');
h.EdgeColor = 'none';
h.FaceColor = 'k';
h.FaceAlpha = 1;
xlabel('Peak voltage [log(uV)]')
ylabel('Probability [%]')
set(gca, 'TickLength', [0 0])
box off
title('Amplitude Distribution')

% max frequency and amplitude vs. time
subplot(3, 4, 7)
scatter(true_pos / ied.fs / 60, peak_val, 2, 'b', 'filled');
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

% ACH 5s
subplot(3, 4, 11)
IED.utils.plotCCG('ccg', ccg, 't', tccg / 1000, 'basepath', basepath,...
    'saveFig', false, 'c', {'k'});
xlabel('Time [s]')
title('Autocorrelogram')

% ACH 10s
subplot(3, 4, 12)
IED.utils.plotCCG('ccg', ccg2, 't', tccg2 / 1000, 'basepath', basepath,...
    'saveFig', false, 'c', {'k'});
xlabel('Time [s]')
title('Autocorrelogram')

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

function markthr(ied,ax)
    % simply add threshold marking to requested ax
    
    % add threshold markers
    if ismember(ied.thrDir,["positive","both"])
        yline(ax, ied.thr(2), '--r');
    end
    if ismember(ied.thrDir,["negative","both"])
        yline(ax, -ied.thr(2), '--r');
    end
end

% EOF