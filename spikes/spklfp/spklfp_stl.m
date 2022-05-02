function stl = spklfp_stl(varargin)

% creates a spike triggered average of lfp by snipping segements of lfp
% sorrounding a subset of spikes. the spikes are selected such that they
% are maximally spaced.

% INPUT
%   basepath    char. fullpath to recording folder {pwd}
%   sig         lfp signal
%   spktimes    cell array of spike times
%   fs          sampling frequency of lfp data
%   mapWin      2 x 1 numeric describing the segment length of lfp signal 
%               to snip arround each spike [s]
%   saveVar     logical  / char. if char then the stl file will be named as
%               the char with the prefix basename
%   graphics    logical {true}
%
% CALLS
%   Sync
%   SyncMap
%
% TO DO LIST
%   snipFromBinary to replace Sync?
%
% 01 apr 22 LH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addParameter(p, 'basepath', pwd, @ischar)
addParameter(p, 'sig', [], @isnumeric)
addParameter(p, 'spktimes', {}, @iscell)
addParameter(p, 'fs', 1250, @isnumeric)
addParameter(p, 'mapWin', [-0.5 0.5], @isnumeric)
addParameter(p, 'graphics', true, @islogical)
addParameter(p, 'saveVar', true)

parse(p, varargin{:})
basepath        = p.Results.basepath;
sig             = p.Results.sig;
spktimes        = p.Results.spktimes;
fs              = p.Results.fs;
mapWin          = p.Results.mapWin;
graphics        = p.Results.graphics;
saveVar         = p.Results.saveVar;

if isempty(mapWin)
    mapWin = [-0.5, 0.5];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% file params
[~, basename] = fileparts(basepath);
if saveVar
    if islogical(saveVar)
    saveVar = 'stl';        
    end
    stlname = [basename, '.', saveVar];
end

% snipping params
winLength = diff(mapWin);
nbinsMap = floor(fs * winLength / 2) * 2 + 1;    % must be odd
maxnspks = 250;
minnspks = 30;

% fft params
fmin = max([0.5, 1 / winLength]);
faxis = [fmin : 0.2 : 100];
faxis = logspace(0, 2, 200);
fftWin = hann(2 ^ (nextpow2(fs) - 1));
noverlap = floor(0.125 * fs);

% sig prep
sig = sig(:);
recDur = length(sig) / fs;                          % [s]
tstamps = [1 / fs : 1 / fs : recDur]';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% snip lfp segments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% loop through cells
for iunit = 1 : length(spktimes)

    bz_Counter(iunit, length(spktimes), 'spk-triggered LFP')

    % select spikes
    nspks = min(maxnspks, length(spktimes{iunit}));
    if nspks < minnspks
        stl.lfp{iunit} = nan(maxnspks, nbinsMap);
        continue
    end
    spkidx = floor(linspace(1, length(spktimes{iunit}), maxnspks));

    % extract win of lfp sorrounding each cell
    [r, i] = Sync([tstamps sig], spktimes{iunit}(spkidx), 'durations', mapWin);
    stl.lfp{iunit} = SyncMap(r, i, 'durations', mapWin,...
        'nbins', nbinsMap, 'smooth', 0);
end

% cat to 3d mat of spike x timebin x unit
stl.lfp = cat(3, stl.lfp{:});

% calc fft per cell
stlAvg = squeeze(mean(stl.lfp, 1, 'omitnan'));
goodidx = all(~isnan(stlAvg));
stl.psd = nan(length(faxis), length(spktimes));
stl.psd(:, goodidx) = pwelch(stlAvg(:, goodidx), fftWin, noverlap, faxis, fs);

% organize and save
stl.info.runtime = datetime(now, 'ConvertFrom', 'datenum');
stl.info.maxnspks = maxnspks;
stl.info.minnspks = minnspks;
stl.faxis = faxis;

if saveVar
    save(fullfile(basepath, stlname), 'stl', '-v7.3')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if graphics
   
    setMatlabGraphics(false)
    fh = figure;
    th = tiledlayout(2, 2, 'TileSpacing', 'Compact');
    title(th, stlname, 'Interpreter', 'none')    
    
    % spike triggered lfp
    nexttile
    xval = linspace(mapWin(1), mapWin(2), nbinsMap);
    ph = plot(xval, mean(stlAvg, 2, 'omitnan'));
    axis tight
    hold on
    plot([0 0], ylim, '--k')
    ylabel('LFP [mV]')
    xlabel('Time [s]')
    title('Mean STL')

    % cell map of psd
    nexttile(2)
    imagesc(faxis, 1 : length(spktimes), stl.psd')
    set(gca, 'xscale', 'log')
    xlim([faxis(1), faxis(end)])
    xlabel('Frequency [Hz]')
    ylabel('Cell No.')
    title('STL PSD per Cell')

    % spike triggered norm psd averaged across cells
    nexttile
    yval = stl.psd ./ sum(stl.psd, 1, 'omitnan');
    ph = plot(faxis, mean(yval, 2, 'omitnan'), 'LineWidth', 3);
    set(gca, 'xscale', 'log')
    xlim([faxis(1), faxis(end)])
    xlabel('Frequency [Hz]')
    ylabel('PSD')
    title('Mean STL PSD')
    
    % psd of stl from all spikes regardless of cell
    nexttile
    yval = squeeze(mean(mean(stl.lfp, 1, 'omitnan'), 3, 'omitnan'));
    [pow, ~] = pwelch(yval, fftWin, noverlap, faxis, fs);
    ph = plot(faxis, pow, 'LineWidth', 3);
    set(gca, 'xscale', 'log')
    xlim([faxis(1), faxis(end)])
    xlabel('Frequency [Hz]')
    ylabel('PSD')
    title('PSD of all Spikes')

   % save figure
    figpath = fullfile(basepath, 'graphics', 'spklfp');
    mkdir(figpath)
    figname = fullfile(figpath, stlname);
    export_fig(figname, '-tif', '-transparent', '-r300')
end


end

% EOF
