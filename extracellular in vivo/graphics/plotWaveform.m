function plotWaveform(avgwv, stdwv, c)

% plots mean +- std waveform of all channels
%
% INPUT
%   avgwv       average waveform [nchans x nsamps] 
%   stdwv       std of waveform [nchams x nsamps]
%   c           color of plot
%
% SEE ALSO
%   plotCluster.m
% 
% 04 dec 18 LH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nargs = nargin;
if nargs < 1 || isempty(avgwv)
    error('average waveform not specified')
end
if nargs < 2 || isempty(stdwv)
    warning('std of waveform not specificed. plotting average only')
    stdwv = [];
end
if nargs < 3 || isempty(c)
    c = 'k';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[nchans, nsamps] = size(avgwv);

for j = 1 : nchans
    offset = j * 100;
    if ~isempty(stdwv)
        errbounds = [avgwv(j, :) + stdwv(j, :);...
            avgwv(j, :) - stdwv(j, :)];
        p = patch([1 : nsamps, nsamps : -1 : 1],...
            [errbounds(1, :), errbounds(2, end : -1 : 1)], c);
        p.EdgeColor = 'none';
        p.FaceAlpha = 0.5;
        set(p, 'YData', get(p, 'YData') - offset);
    end
    hold on
    l = plot(avgwv(j, :), 'lineWidth', 1, 'Color', c);
    set(l, 'YData', get(l, 'YData') - offset);
end
axis off

end

% EOF


