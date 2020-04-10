function [linet] = lineDetect(varargin)

% wrapper for firing rate functions
% 
% INPUT
%   x           signal 
%   fs          sampling frequency
%   bandpass    low and high cutoffs {[50 70]}
%   mode        filtering mode, {'iir'}
%   phase0      the phase to detect {-pi/2}; zero phase corresponds to peaks (do
%               not use for subtraction - temporal variability is then reflected in
%               spikes at cycle peaks)
%   graphics    for a short segment

% OUTPUT
% 	line        zero-crossing times [samples]      
%
% CALLS
%   firfilt
%   separators
% 
% 29 dec 19 LH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addOptional(p, 'x', []);
addOptional(p, 'fs', 1250);
addOptional(p, 'bandpass', [50 70]);
addOptional(p, 'mode', 'iir', @ischar);
addOptional(p, 'phase0', -pi/2);
addOptional(p, 'graphics', true, @islogical);

parse(p, varargin{:})
x = p.Results.x;
fs = p.Results.fs;
bandpass = p.Results.bandpass;
mode = p.Results.mode;
phase0 = p.Results.phase0;
graphics = p.Results.graphics;

x = double(x(:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filtered = filterLFP(x, 'fs', fs, 'passband', bandpass, 'order', 6,...
    'type', 'butter', 'dataOnly', true, 'graphics', false,...
    'saveVar', false);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  neg >> pos zero-crossings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

phs = angle(hilbert(filtered) ); % zero = pos. peak
linet = find([0; diff( phs > phase0)] == 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if graphics
    
    tidx = (1 : round(fs / 10));    % 100 ms
    idx = round(fs / 2) + tidx;     % starting from 500 ms
    t = idx / fs;
    lidx = find(linet > idx(1) & linet < idx(end));
    
    figure
    subplot(1, 2, 1)
    plot(t, x(idx), 'b')
    hold on
    plot(t, filtered(idx), 'r')
    legend('raw', 'filtered')
    line([linet(lidx) linet(lidx)] / fs, ylim, 'Color', 'k', 'HandleVisibility','off')
    box off
    axis tight
    set(gca, 'TickLength', [0 0])
    xlabel('Time [s]')
    
    subplot(1, 2, 2)
    cycledur = diff(linet / fs * 1000);
    h = histogram(diff(linet / fs * 1000));
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

