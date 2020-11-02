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
%               selected grp {1}.
%   vis         char. figure visible {'on'} or not ('off')
%
% CALLS
%   none
%
% OUTPUT
%   fepsp       struct with fields described below
%
% TO DO LIST
%   # Lior Da Marcas take over
%
% 16 oct 20 LH   UPDATES

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

parse(p, varargin{:})
fepsp = p.Results.fepsp;
basepath = p.Results.basepath;
dt = p.Results.dt;
force = p.Results.force;
saveVar = p.Results.saveVar;
saveFig = p.Results.saveFig;
graphics = p.Results.graphics;
vis = p.Results.vis;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get params from fepsp struct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~, basename] = fileparts(basepath);
fepspname = [basename '.fepsp.mat'];
    
% try to load file if not in input
if isempty(fepsp)
    if exist(fepspname, 'file')
        load(fepspname)
    end
end
if isfield('fepsp', 'amp') && ~force
    load(fepspname)
    return
end

fs = fepsp.info.fs;
tstamps = fepsp.tstamps;
spkgrp = fepsp.info.spkgrp;
nspkgrp = length(spkgrp);
nfiles = length(fepsp.intens);
protocol = fepsp.info.protocol;

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
        switch fepsp.info.recSystem
            case 'oe'
                % correct stim frequency
                ts = fepsp.info.stimTs;
                ts = mean(ts(ts < 500)) / fs * 1000;
            case 'wcp'
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
% ampcell       2d array (tetrode x intensity) where each cell contains the
%               amplitude/s for each trace
fepsp.ampcell   = cell(nspkgrp, nfiles);
% amp           2d (io) or 3d (stp) mat (tetrode x intensity x stim) of
%               amplitude averaged across traces
fepsp.amp   = nan(nspkgrp, nfiles, nstim);
% ampNorm       2d array of normalized amplitudes. for each trace the
%               responses are normalized to first response. these
%               normalized amplitudes. for
%               stp only
fepsp.ampNorm   = nan(nspkgrp, nfiles, nstim);
% facilitation  2d mat of average maximum normalized response. for stp only
fepsp.facilitation = nan(nspkgrp, nfiles);

fepsp = rmfield(fepsp, 'ampNorm');

for j = 1 : nspkgrp
    for i = 1 : nfiles
        fepsp.traceAvg(j, i, :) = mean(fepsp.traces{j, i}, 2);
        switch protocol
            case 'io'
%                 Generated bad results
%                 if abs(min(fepsp.traces{j, i}(:))) > max(fepsp.traces{j, i}(:))
%                     % I took this method from getLFP line 178, however
%                     % changed the direction as it seemed weird the
%                     % opposite. However, the syncronus/ asyncronus/
%                     % inhibitory responce might affect this, and maybe with
%                     % fepsp in PYR layer we would expect opposite direction.
%                     % this will affect calculation because of the diff sensetivity of
%                     % 1st local max and 1st local min, so we need to invert.
%                     fprintf('Inverting traces in cell (%d,%d)\n',j,i)
%                     fepsp.traces{j, i} = -fepsp.traces{j, i};
%                 end
                fepsp.waves{j, i} = fepsp.traces{j, i}(wvwin(1) : wvwin(2), :);
                % find first local minima & maxima in window. since we
                % expect the minima to be relativly big, ignore any minima
                % with small Prominence (<0.3). if the responce is smaller
                % then that, it is practicly 0.
                [ValueCheckMax,PerTrace1Max] = max(islocalmax(fepsp.waves{j, i},1),[],1);
                [ValueCheckMin,PerTrace1Min] = max(islocalmin(fepsp.waves{j, i},1,'MinProminence',0.3),[],1);
                PerTrace1Max(ValueCheckMax == 0) = 1; 
                %if no local maxima with any Prominence, take the 1st sample in window
                PerTrace1Min(ValueCheckMin == 0) = round((fs/1000)*10-2*(fs/1000));
                %if no local minima with minimal Prominence, take the sample corresponding to 10ms (expected responce pos)
                PerTrace1MinIND = sub2ind(size(fepsp.waves{j, i}),PerTrace1Min,1:size(fepsp.waves{j, i},2));
                PerTrace1MaxIND = sub2ind(size(fepsp.waves{j, i}),PerTrace1Max,1:size(fepsp.waves{j, i},2));
                fepsp.ampcell{j, i} = abs(fepsp.waves{j, i}(PerTrace1MaxIND) - fepsp.waves{j, i}(PerTrace1MinIND));
                fepsp.wavesAvg(j, i, :) = mean(fepsp.waves{j, i}, 2);
                %fepsp.ampcell{j, i} = range(fepsp.waves{j, i});
                %fepsp.amp(j, i) = mean(fepsp.ampcell{j, i});
                
                %Redo on mean wave
                TheMin{j,i} = find(islocalmin(squeeze(fepsp.wavesAvg(j, i, :)),'MinProminence',0.3),1);
                TheMax{j,i} = find(islocalmax(squeeze(fepsp.wavesAvg(j, i, :))),1);
                if isempty(TheMin{j,i})
                    TheMin{j,i} = round((fs/1000)*10-2*(fs/1000)); 
                    %if no local minima with minimal Prominence, take the sample corresponding to 10ms (expected responce pos)
                end
                if isempty(TheMax{j,i})
                    TheMax{j,i} = 1; %if no local maxima with any Prominence, take the 1st sample in window
                end
                fepsp.amp(j, i) = abs(fepsp.wavesAvg(j, i, TheMax{j,i}) - fepsp.wavesAvg(j, i, TheMin{j,i}));
                
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
                fh = figure('Visible', vis);
                suptitle(sprintf('T#%d', i))
                subplot(1, 2, 1)
                plot(fepsp.tstamps(wvwin(1) : wvwin(2)), squeeze(fepsp.wavesAvg(i, :, :))')
                hold on
                SmallTSwindow = fepsp.tstamps(wvwin(1) : wvwin(2));
                wavesAvgPos = sub2ind(size(fepsp.wavesAvg),ones(1,length([TheMin{i,:} TheMax{i,:}]))*i,repmat(1:size(fepsp.wavesAvg,2),1,2),[TheMin{i,:} TheMax{i,:}]);
                plot(SmallTSwindow([TheMin{i,:}]),squeeze(fepsp.wavesAvg(wavesAvgPos(1:(end/2)))'),'*')
                plot(SmallTSwindow([TheMax{i,:}]),squeeze(fepsp.wavesAvg(wavesAvgPos((end/2+1):end))'),'O')
                axis tight
                yLimit = [min([fepsp.wavesAvg(:)]) max([fepsp.wavesAvg(:)])];
                ylim(yLimit)
                xlabel('Time [ms]')
                ylabel('Voltage [mV]')
                legend([split(num2str(sort(fepsp.intens)));{'Amp measure point 1_s_t min'};{'Amp measure point 1_s_t max'}])
                box off
                
                subplot(1, 2, 2)
                %ampmat = cell2nanmat(fepsp.ampcell(i, :));
                %boxplot(ampmat, 'PlotStyle', 'traditional')
                bar(fepsp.amp(i,:))
                ylim([min(horzcat(fepsp.ampcell{:})) max(horzcat(fepsp.ampcell{:}))])
                xticklabels(split(num2str(sort(fepsp.intens))))
                xlabel('Intensity [uA]')
                ylabel('Amplidute [mV]')
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
            figname = [figpath '\fepsp_t' num2str(i)];
            export_fig(figname, '-tif', '-r300', '-transparent')
        end
    end
end

end

% EOF