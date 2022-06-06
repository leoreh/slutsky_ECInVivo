function [revIdx] = tsa_detectRev(varargin)

% detects revolutions (phase crossings in a frequency band)
% 
% INPUT
%   sig         numeric. signal 
%   fs          numeric. sampling frequency
%   bandpass    numeric. low and high cutoffs {[50 70]}
%   p0          numeric. the phase to detect in radians. {-pi/2}.
%   graphics    logical. for a short segment

% OUTPUT
% 	revIdx      numeric. phase-crossing times [samples]      
%
% CALLS
%   firfilt
% 
% jul 18 LH & ES

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addOptional(p, 'sig', [], @isnumeric);
addOptional(p, 'fs', 1250, @isnumeric);
addOptional(p, 'bandpass', [45 55], @isnumeric);
addOptional(p, 'p0', -pi/2, @isnumeric);
addOptional(p, 'graphics', true, @islogical);

parse(p, varargin{:})
sig = p.Results.sig;
fs = p.Results.fs;
bandpass = p.Results.bandpass;
p0 = p.Results.p0;
graphics = p.Results.graphics;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% work
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sig = double(sig(:));

% filter
filtered = filterLFP(sig, 'fs', fs, 'passband', bandpass, 'order', 6,...
    'type', 'butter', 'dataOnly', true, 'graphics', false,...
    'saveVar', false);


%  neg >> pos zero-crossings
phs = angle(hilbert(filtered)); % zero = pos. peak
revIdx = find([0; diff(phs > p0)] == 1);
revIdx = revIdx(:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if graphics
    
    tidx = (1 : round(fs / 10));    % 100 ms
    idx = round(fs / 2) + tidx;     % starting from 500 ms
    t = idx / fs;
    lidx = find(revIdx > idx(1) & revIdx < idx(end));
    
    figure
    subplot(1, 2, 1)
    plot(t, sig(idx), 'b')
    hold on
    plot(t, filtered(idx), 'r')
    legend('raw', 'filtered')
    line([revIdx(lidx) revIdx(lidx)] / fs, ylim, 'Color', 'k', 'HandleVisibility','off')
    box off
    axis tight
    set(gca, 'TickLength', [0 0])
    xlabel('Time [s]')
    
    subplot(1, 2, 2)
    h = histogram(diff(revIdx / fs * 1000));
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

