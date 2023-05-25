function plot_wv(varargin)

% plots mean +- std waveform of all channels horizontal (concatenated) or
% vertical
%
% INPUT
%   wv          average waveform [nchans x nsamps] 
%   wv_std      std of waveform [nchams x nsamps]
%   clr         color of plot
%   orient      channels horizontal {horz} or vertical (vert)
%   sbar        plot scale bar {1} or not (0)
%   fs          sampling rate
%
% 04 dec 18 LH  updates:
% 22 jan 19     horizontal plot
% 08 may 19     offset according to trace
% 06 oct 19     offset according to known value
% 18 jun 22     clean up

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addOptional(p, 'wv', []);
addOptional(p, 'wv_std', []);
addOptional(p, 'clr', 'k', @ischar);
addOptional(p, 'orient', 'horz', @ischar);
addOptional(p, 'fs', 24414.06, @isnumeric);
addOptional(p, 'sbar', true, @islogical);

parse(p,varargin{:})
wv          = p.Results.wv;
wv_std      = p.Results.wv_std;
clr         = p.Results.clr;
orient      = p.Results.orient;
fs          = p.Results.fs;
sbar        = p.Results.sbar;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[nchans, nsamps] = size(wv);
j = 1;
if strcmp(orient, 'horz')
        wv = [wv, NaN(nchans, nsamps / 4)];
        wv = reshape(wv', [], 1);
        if ~isempty(wv_std)
            wv_std = [wv_std, zeros(nchans, nsamps / 4)];
            wv_std = reshape(wv_std', [], 1);
            errbounds = [wv + wv_std, wv - wv_std];
            errbounds(isnan(errbounds)) = 0;
            p = patch([1 : length(errbounds), length(errbounds) : -1 : 1],...
                [errbounds(:, 1); errbounds(end : -1 : 1, 2)]', clr);
            p.EdgeColor = 'none';
            p.FaceAlpha = 0.3;
        end
        hold on
        l = plot(wv, 'lineWidth', 2, 'Color', clr);
%         ylim([-200 200])
else
    for j = 1 : nchans
        % offset = j * (max(wv(j, :)) - min(wv(j, :)));
        offset = 500 * (j - 1);
        if ~isempty(wv_std)
            errbounds = [wv(j, :) + wv_std(j, :);...
                wv(j, :) - wv_std(j, :)];
            p = patch([1 : nsamps, nsamps : -1 : 1],...
                [errbounds(1, :), errbounds(2, end : -1 : 1)], clr);
            p.EdgeColor = 'none';
            p.FaceAlpha = 0.3;
            set(p, 'YData', get(p, 'YData') - offset);
        end
        hold on
        l = plot(wv(j, :), 'lineWidth', 1, 'Color', clr);
        set(l, 'YData', get(l, 'YData') - offset);
    end
end
axis off

if sbar
    % scale bar:
    % time (x-axis), sampling frequency 24414
    % voltage (y-axis), the PZ5 by TDT has an input range of +/- 500 mV,
    % digitized @ 28 bit resolution. Accordingly:
    % signal [V] = signal [bits] * IO range [V] / ADC resolution [nbits]
    % nbits = 2 ^ 27 -1;          % resolution [bits]
    % maxv = 1;                 % maximum voltage [V]
    % bit = maxv / nbits;         % one bit [V]
    % sigv = double(x) * bit * 10 ^ 6
    
    line([0 0], [0 -100], 'Color', 'k')
    line([0 fs / 10 ^ 3], [0 0], 'Color', 'k')
end

end

% EOF


