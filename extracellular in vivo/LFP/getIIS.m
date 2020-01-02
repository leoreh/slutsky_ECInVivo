function [iis] = getIIS(varargin)

% detects inter-ictal spikes from LFP.
% 
% INPUT
%   sig         signal for detection
%   fs          sampling frequency
%   binsize     scalar {60}. for rate calculation
%   marg        scalar {0.1} in [s]. time margin for clipping spikes 
%   thr         scalar {15} in [z scores]. used for detection 
%   basepath    recording session path {pwd}
%   basename    string. if empty extracted from basepath
%   graphics    logical {true}. plot figure
%   saveVar     logical {true}. save variable
%   saveFig     logical {true}. save figure
%   forceA      logical {false}. force analysis even if .mat exists
%
% OUTPUT
%   iis         struct with fields:
%       stamps      n x 2 mat where n is the number of events and the
%                   columns represent the start and end of each event [samps]
%       peaks       the timestamp for each peak [samps]
%       peakPower   the voltage at each peak [uV]
%       <params>    as in input
%
% TO DO LIST
%       # investigate thr
%
% CALLS
%       calcFR     
%
% 02 jan 20 LH.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addParameter(p, 'sig', [], @isnumeric)
addParameter(p, 'fs', 1250, @isnumeric)
addParameter(p, 'binsize', 600, @isnumeric)
addParameter(p, 'marg', [], @isnumeric)
addParameter(p, 'thr', 15, @isnumeric)
addParameter(p, 'basepath', pwd, @isstr);
addParameter(p, 'basename', [], @isstr);
addParameter(p, 'graphics', true, @islogical)
addParameter(p, 'saveVar', true, @islogical);
addParameter(p, 'saveFig', true, @islogical);
addParameter(p, 'forceA', false, @islogical);

parse(p, varargin{:})
sig = p.Results.sig;
fs = p.Results.fs;
marg = p.Results.marg;
thr = p.Results.thr;
binsize = p.Results.binsize;
basepath = p.Results.basepath;
basename = p.Results.basename;
graphics = p.Results.graphics;
saveVar = p.Results.saveVar;
saveFig = p.Results.saveFig;
forceA = p.Results.forceA;

if isempty(marg)
    marg = 0.1 * fs;
end
tstamps = [1 : length(sig)] / fs;

% params
thr = 15;               % for detection [z scores]
marg = 0.1 * fs;        % time margin for clipping spikes

% initialize output
iis = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check if mat already exists
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(basename)
    [~, basename] = fileparts(basepath);
end
filename = [basepath, '\', basename, '.iis.mat'];
if ~forceA
    if exist(filename)
        load(filename)
        fprintf('\n loading %s \n\n', filename)
        return
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prep signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = zscore((sig / rms(sig)) .^ 2);

% investigate threshold
% j = 1;
% idx = 5 : 0.5 : floor(max(x));
% nevents = zeros(1, length(idx));
% for i = 1 : length(idx)
%     nevents(i) = sum(x > idx(i));
% end
% figure
% plot(idx, cumsum(nevents))
% yyaxis right
% plot(idx, cumsum(log10(nevents)))
% axis tight
% 
% histogram((nevents), 50, 'Normalization', 'cdf')

        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% detect signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iie = find(diff(x > thr) > 0);

% select local maximum and clip spike
peak = zeros(length(iie), 1);
pos = zeros(length(iie), 1);
seg = zeros(length(iie), marg * 2 + 1);
for i = 1 : length(iie)
    seg(i, :) = sig(iie(i) - marg : iie(i) + marg);
    [peak(i), pos(i)] = max(abs(seg(i, :)));
    pos(i) = iie(i) - marg + pos(i);
    seg(i, :) = sig(pos(i) - marg : pos(i) + marg);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% power spectrum via wavelet
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fseg = linspace(1, 250, 250);
% pwelch(seg(1, :), [], [], fseg, 1250);
%
% cwt(seg(1, :), 1250)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rate 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[iis.rate, iis.edges, iis.cents] = calcFR(pos, 'winCalc', [1, length(sig)],...
    'binsize', 60, 'smet', 'none', 'c2r', false);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arrange struct output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% bs.stamps = temp;
% bs.peakPos = peakPos;
% bs.peakPower = peakPower;
% bs.dur = dur;
% bs.iinterval = iinterval;
% bs.fs = fs;
% bs.binsize = binsize;
% bs.vars = vars;
% bs.lRat = lRat;
% bs.iDist = iDist;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if saveVar
    save(filename, 'bs')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if graphics
    fh = figure;
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    
    % raw and iis
    subplot(3, 3, 1 : 2)
    plot(tstamps / 60, sig)
    yyaxis right
    plot(tstamps / 60, x, 'k')
    hold on
    plot(xlim, [thr thr], '--y')
    plot([pos pos] / fs / 60, [-10 -1], '--g', 'LineWidth', 2)
    xlim([120 130])
    
    % iis rate
    subplot(3, 3, 4 : 5)
    plot(iis.cents / fs / 60, iis.rate, 'k', 'LineWidth', 1)
    
    % iis waveforms
    subplot(3, 3, 3)
    xstamps = [1 : size(seg, 2)] / fs;
    plot(xstamps, seg')
    hold on
    axis tight
    stdshade(seg, 0.5, 'k', xstamps)
    
    if saveFig
        figname = [basename '_BS'];
        export_fig(figname, '-tif', '-transparent')
        % savePdf(figname, basepath, ff)
    end
end
end

% EOF

