function [sig, tsaSig, revIdx] = tsa_filter(varargin)

% filters a signal in the time domain by time synchroneous averaging of a
% periodic interference. finds the length of each revolution (period) by
% filtering the signal in the relevant passband (see tsa_detectRev) and
% averages the signal in these periods to isolate the interfernece waveform
% (tsa). the tsa is than removed from each period of the signal. 

% INPUT
%   sig         numeric. signal 
%   revIdx      numeric. phase-crossing times [samples]  
%   ma          logical. use moving average (true) or weighted average
%               false)
%   tw          logical. use time warpping {true} or clipping (false) 
%               when removing the tsa
%   fs          numeric. sampling frequency
%   bandpass    numeric. low and high cutoffs {[45 55]}
%   p0          numeric. the phase to detect in radians. {-pi/2}.
%   graphics    logical.

% OUTPUT
% 	sig         numeric. signal without the piece-wise revIdxar effect of the
%               interference
%   tsaSig      numeric. tsa (average interference) of the signal
%
% CALLS
%   tsa_detectRev
% 
% jul 18 LH & ES
%
% algorithmic notes:
%
% (1) choice of update method:
%   to track variations in the noise (the interference), the MA method is
%   optimal. However, it also tracks variations in the signal. For instance,
%   if there is a MA of 60 seconds and there is a short deflection (e.g.
%   pulse) of 60 ms, the contamination factor is 1:1000. But if there is a
%   pulse train with 50% duty cycle, the contamination will be 50%. In the
%   case of a zero-base line signal with a pulse train of amplitude V, the
%   mean of the xlta will thus be V/2 (instead of close to zero in the
%   weighted average method), which will induce a huge shift in the signal.
%
% (2) signal mean:
%   in case of a zero-mean raw signal, the w-method will work perfectly even in
%   the presence of temporary deflections. However, in a DC-coupled signal
%   with a mean M, the mean of the xlta will also be M. This will result in
%   the revIdx-removed waveform having a zero-mean. If the raw signal has
%   zero-baserevIdx but occasional deflections of amplitude N (e.g. the mean
%   is 0.5 and the max value is 1, i.e. a duty cycle d of 50%), then
%   the mean of the xlta will be 0.5. Subtracting will yield a zero-mean
%   signal, with a baserevIdx of -0.5 and a max value of 0.5.
%   this problem is solved by adding the DC to the xlta. Thus, the formula
%   is
%
%       y(t) = x(t) - n(t) + myu(x)
%
%   in the above example, adding the mean of 0.5 will shift the signal back
%   to have zero baserevIdx and max of 1.
%
%   the situation is less than ideal when the duty cycle is small and the
%   method is MA. For instance, assume dc is 1%, the signal baserevIdx M is
%   e.g. -1, and the max deflection N is e.g. 2. Then, the
%   signal mean is
%           M * ( 1 - dc ) + ( M + N ) * dc =
%           M + N * dc =
%           -0.98
%   if the instantaneous xlta is for a period in which the mean is M, then
%   the error will be
%
%           -M + [ M + N * dc ] = N * dc
%
%   which, in our example, 2*0.01 = 0.02 positive shift. Of course, this is
%   still much better than not removing the mean
%
% (3) the I.C. should be computed for a period without any major deflections.
%
% (4) rarely-sampled points:
%   are treated by a weighted average between the sample, the last
%   fully-sampled point, and the cyclically next fully-sampled point. note
%   that this is relevant only for the w-average mode

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addOptional(p, 'sig', [], @isnumeric);
addOptional(p, 'revIdx', [], @isnumeric);
addOptional(p, 'ma', true, @islogical);
addOptional(p, 'tw', true, @islogical);
addOptional(p, 'fs', 1250, @isnumeric);
addOptional(p, 'bandpass', [45 55], @isnumeric);
addOptional(p, 'p0', -pi/2, @isnumeric);
addOptional(p, 'graphics', true, @islogical);

parse(p, varargin{:})
sig = p.Results.sig;
revIdx = p.Results.revIdx;
ma = p.Results.ma;
tw = p.Results.tw;
fs = p.Results.fs;
bandpass = p.Results.bandpass;
p0 = p.Results.p0;
graphics = p.Results.graphics;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% gen params
if graphics
    orig = sig(:);  % for graphics and debugging
end
sig = sig(:);
dc = mean(sig);
w = 0.9998;
nsamps = length(sig);

% revolution indices
if isempty(revIdx)
    revIdx = tsa_detectRev('sig', sig, 'fs', fs, 'graphics', false,...
        'bandpass', bandpass, 'p0', p0);
end
dRev = diff(revIdx);
avgRev = round(mean(dRev));
maxRev = max(dRev);
n = ceil(1 / (1 - w));
win = [0, min([max(revIdx) - maxRev, maxRev * n])];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calc tsa
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize tsa mat and weight vectors
if tw
    xmat = zeros(n, avgRev);
    wvec = ones(1, avgRev) * n;
else
    xmat = zeros(n, maxRev);

    % handle rarely-sampled points by calculating weighted mean
    dlic = dRev(1 : n);
    edges = 0.5 : 1 : max(dlic) + 0.5;
    h = histcounts(dlic, edges);
    h(end) = [];
    Fx = (sum(h) - cumsum(h));
    wvec = [n, Fx];
    wvec(end : maxRev) = zeros(1, maxRev - length(Fx));
    
end

% prepare initial tsa
for irev = 1 : win(2) / maxRev
    t1 = revIdx(irev);
    t2 = revIdx(irev + 1) - 1;
    if t2 > nsamps
        t2 = nsamps;
    end
    xmat(irev, :) = tsa_prepseg(sig(t1 : t2), size(xmat, 2), tw);
end
% tsaSig = sum(xmat, 1, 'omitnan') ./ wvec;
% tsaSig(isnan(tsaSig)) = 0;
tsaSig = mean(xmat);
[n, m] = size(xmat);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% adaptive subtraction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for irev = 1 : length(revIdx) - 1
    t1 = revIdx(irev);
    t2 = revIdx(irev + 1) - 1;
    nx = t2 - t1 + 1;
    if nx > maxRev
        t2 = maxRev;
        sig(t1 : t2)  = sig(t1 : t2) - tsaSig(1 : length(t2))';
        continue
    end
    if tw
        idx = 1 : (avgRev - 1) / (nx - 1) : avgRev;
        tsax = tsaSig(round(idx));
    else
        idx = 1 : nx;
        if nx > length(tsaSig)
            tsax = [tsaSig, zeros(1, nx - length(tsaSig))];
        else
            tsax = tsaSig(idx);
        end
    end
    sig(t1 : t2) = sig(t1 : t2) - tsax';

    % update average 
    if irev > win(2) / maxRev
        x = tsa_prepseg(sig(t1 : t2), m, tw);
        if length(x) > length(tsaSig)
            x = x(1 : length(tsaSig));
        end
        if ma
            ihat = mod(irev - 1, n) + 1;
            tsaSig = tsaSig - (1 - w) * xmat(ihat, :) + (1 - w) * x;
            xmat(ihat, :) = x;
            if ~tw
                wvec = sum(xmat ~= 0)';
            end
        else
            if tw
                tsaSig = w * tsaSig + (1 - w) * x;
            else
                tsaSig(idx) = w * tsaSig(idx) + (1 - w) * x(idx);
                wvec(idx) = wvec(idx) + ones(1, nx);
            end
        end
    end
    if t2 >= nsamps
        break
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if graphics
    setMatlabGraphics(false)
    fh = figure;
    
    % example of the clean and raw signal
    subplot(2, 2, 1)
    tstart = round(length(sig) / 2);
    nRev = 5;
    sigIdx = [tstart : nRev * length(tsaSig) + tstart];
    tstamps = sigIdx / fs;
    lidx = find(revIdx > sigIdx(1) & revIdx < sigIdx(end));
    plot(tstamps, orig(sigIdx), 'k')
    hold on
    plot(tstamps, sig(sigIdx), 'b')
    line([revIdx(lidx) revIdx(lidx)] / fs, ylim, 'Color', 'k', 'HandleVisibility','off')
    legend('raw', 'filtered')
    box off
    axis tight
    set(gca, 'TickLength', [0 0])
    xlabel('Time [s]')
    ylabel('Amplitude')
    
    % pwelch on 30 min of signal
    sigLim = min([30 * 60 * fs, nsamps]);
    sigIdx = [1 : sigLim];    
    subplot(2, 2, 2)
    fftwin = hann(2 ^ (nextpow2(2 * fs) - 1));
    noverlap = floor(0.25 * fs);
    faxis = [0.2 : 0.2 : 220];
    [pow_orig, ~] = pwelch(orig(sigIdx), fftwin, noverlap, faxis, fs);
    [pow_clean, ~] = pwelch(sig(sigIdx), fftwin, noverlap, faxis, fs);
    ph = plot(faxis, pow_orig ./ sum(pow_orig, 2), 'k', 'LineWidth', 2);
    hold on
    ph = plot(faxis, pow_clean ./ sum(pow_clean, 2), 'b', 'LineWidth', 2);
    legend('raw', 'filtered')
    xlabel('Frequency [Hz]')
    ylabel('Norm PSD')
    
    % tsa
    subplot(2, 2, 3)
    plot(tsaSig)
    xlabel('sample')
    ylabel('amplitude')
    title('average interference')
    axis tight
    
    % histogram of cycle duration
    subplot(2, 2, 4)
    h = histogram(diff(revIdx / fs * 1000), 100,...
        'Normalization', 'probability');
    xlabel('cycle duration [ms]')
    ylabel('counts')
    h.EdgeColor = 'none';
    h.FaceColor = 'k';
    box off
    axis tight
    set(gca, 'TickLength', [0 0])
end

end

% EOF

% -------------------------------------------------------------------------
% free call
% graphics = true;
% ma = true;
% tw = false;
% bandpass = [45 55];
% p0 = -pi/2;
% revIdx = [];
% fs = 1250; 
% idx = [5 * 60 * 60 * fs : 10 * 60 * 60 * fs];
% orig = EMG(idx);
% sig = orig;