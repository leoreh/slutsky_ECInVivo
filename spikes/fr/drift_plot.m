function drft = drift_plot(drft, axh)

% plots PV correlations and drift. see drift_calc
%
% INPUT:
%   drift           struct. see drift_calc
%   axh             axis handle. if empty will create new figure
%
% DEPENDENCIES:
%   none
%
% 22 may 24 LH      based on Lee's code

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
    axh = [];
end

if isempty(axh)
    setMatlabGraphics(true)
    fh = figure;
    set(fh, 'WindowState', 'maximized');
    tlayout = [1, 2];
    th = tiledlayout(tlayout(1), tlayout(2));
    th.TileSpacing = 'tight';
    th.Padding = 'none';
    set(fh, 'DefaultAxesFontSize', 16);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xval = [1 : length(drft.m_corr)];

hold on
plot(xval, drft.dt_corr, '.', 'Color', [.5 .5 .5]);
plot(xval, drft.m_corr, '.k', 'MarkerSize', 15);
plot(xval, drft.lin_coef(2) + drft.lin_coef(1) * xval,'b');

xticks([1 : 3600 / drft.info.winsize : xval(end)]);
xticklabels(1 : xval(end)); 

title(sprintf('%Drift Rate = %.3f', drft.drate),...
    'interpreter', 'none', 'FontSize', 20);
xlabel(sprintf('\x394 Time [h]'))
ylabel('PV correlation');
axis tight



end

% EOF

