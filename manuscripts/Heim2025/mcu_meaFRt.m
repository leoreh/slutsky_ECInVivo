function mcu_meaFRt(v, alt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MEA Firing Rate Time Analysis
% 
% Inputs:
%   v - cell array of data structures (if empty, loads data automatically)
%   alt - analysis option:
%       1: Plot MFR+SD for both groups on single plot
%       2: Plot FR by percentiles
%       3: Plot FR-fit unit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load data if not provided
if isempty(v)
    grps = {'mea_bac'; 'mea_mcuko'};
    grpLbls = {'Control'; 'MCU-KO'};
    vars = {'frr', 'st_brst'};
    
    for iGrp = 1 : 2
        basepaths = mcu_sessions(grps{iGrp});
        v{iGrp} = basepaths2vars('basepath', basepaths, 'vars', vars);
    end
else
    grpLbls = {'Control'; 'MCU-KO'};
end

% Get perturbation onset and create time vector
idxPert = v{1}(1).frr.info.idxPert(1);
binSize = v{1}(1).frr.info.binSize(1);
t = ((1:size(v{1}(1).frr.fr, 2)) - idxPert) * binSize / 60 / 60 * 3;
tPert = t(idxPert);

% Run selected analysis
switch alt
    case 1
        plot_mfr_both_groups(v, t, tPert, grpLbls);
    case 2
        plot_fr_percentiles(v, t, tPert, grpLbls);
    case 3
        plot_fr_fit_unit(v, t, tPert, grpLbls);
    otherwise
        error('Invalid alt option. Use 1, 2, or 3.');
end

end

function plot_mfr_both_groups(v, t, tPert, grpLbls)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot MFR+SD for both groups on single plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set graphics
[hFig, hAx] = plot_set('Firing Rate (% BSL)');
clr = mcu_clr;
hold on;

% Plot each group
for iGrp = 1 : 2
    % Cat data
    frr = catfields([v{iGrp}(:).frr], 1);
    uGood = frr.uGood;
    
    % Calculate MFR
    fr = frr.fr(uGood, :);
    frBsl = frr.mdlF.frBsl(uGood);
    frMat = (fr ./ frBsl) * 100;
    
    % Plot MFR + SD
    hGrp = plot_stdShade('dataMat', frMat,...
        'xVal', t, 'hAx', hAx, ...
        'clr', clr.grp(iGrp, :), 'alpha', 0.5);
    set(hGrp, 'DisplayName', sprintf('%s (n=%d)', grpLbls{iGrp}, sum(uGood)));
end

% Formatting
hLgd = legend;
hLgd.Location = 'southeast';
hAx.FontSize = 16;

% Mark perturbation and baseline
ylim([0 140])
plot_pert(hAx, tPert, clr.bac);
yline(100, '--', '', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1,...
    'HandleVisibility', 'off');

plot_axSize('hFig', hFig, 'szOnly', false, 'axShape', 'wide',...
    'axHeight', 300, 'flgPos', true);

% Save
fname = 'MEA ~ FRtime';
lme_save('hFig', hFig, 'fname', fname, 'frmt', {'svg', 'mat'});

end

function plot_fr_percentiles(v, t, tPert, grpLbls)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot FR by percentiles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Hard-coded parameters
iGrp = 1;
mdlPrfx = 'mdlF';
nPrc = 3;
grpVar = 'fr';
grpVar = 'bspks';

% Cat data
frr = catfields([v{iGrp}(:).frr], 1);
brst = catfields([v{iGrp}(:).brst], 2);
uGood = frr.uGood;

% Get grouping data
if strcmp(grpVar, 'fr')
    grpData = frr.(mdlPrfx).frBsl(uGood);
else
    grpData = brst.bspks(1, uGood);
end

% Create percentiles
alpha = 1.5; 
p = linspace(0, 1, nPrc + 1).^alpha;
percEdges = prctile(grpData, 100 * (1 - p));
percEdges = sort(percEdges);

% Set graphics
[hFig, hAx] = plot_set('Firing Rate (% BSL)');
clr = mcu_clr;
clrPrc = bone(nPrc + 2);
clrPrc = clrPrc(1 : end - 1, :);
hold on;

% Plot each percentile
for iPrc = 1 : nPrc
    % Select units in this percentile
    if iPrc == 1
        uPrc = grpData <= percEdges(iPrc + 1);
    else
        uPrc = grpData > percEdges(iPrc) & grpData <= percEdges(iPrc + 1);
    end
    uPlt = find(uGood);
    uPlt = uPlt(uPrc);

    % Calculate MFR from selected units
    if strcmp(mdlPrfx, 'mdlF')
        fr = frr.fr(uPlt, :);
    else
        fr = frr.frMdl(uPlt, :);
    end
    frBslPerc = frr.(mdlPrfx).frBsl(uPlt);

    % Plot MFR + SD
    frMat = (fr ./ frBslPerc) * 100;
    hPtl = plot_stdShade('dataMat', frMat,...
        'xVal', t, 'hAx', hAx, ...
        'clr', clrPrc(iPrc, :), 'alpha', 0.5);

    % Create legend text with range
    if strcmp(grpVar, 'fr')
        strPrc = sprintf('%.0f-%.0f', percEdges(iPrc), percEdges(iPrc + 1));
        txtLgd = sprintf('%s Hz (n=%d)', strPrc, length(uPlt));
    else
        strPrc = sprintf('%s=%.1f-%.1f', 'P(S\inB)', percEdges(iPrc), percEdges(iPrc + 1));
        txtLgd = sprintf('%s (n=%d)', strPrc, length(uPlt));
    end
    set(hPtl, 'DisplayName', txtLgd);
end

% Formatting
hLgd = legend('Interpreter', 'tex');
if iGrp == 1
    hLgd.Location = 'southeast';
    % ylim([0 200]);
else
    hLgd.Location = 'northeast';
    hLgd.Position(2) = hLgd.Position(2) - 0.05;
    % ylim([0 180])
end
title(sprintf('%s', grpLbls{iGrp}));

% Mark perturbation and baseline
plot_pert(hAx, tPert, clr.bac);
yline(100, '--', '', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1,...
    'HandleVisibility', 'off');
plot_axSize('hFig', hFig, 'szOnly', false, 'axShape', 'wide', 'axHeight', 300);

% Save
% fname = ['MEA ~ FRtime_Prctile_', grpLbls{iGrp}];
% lme_save('hFig', hFig, 'fname', fname, 'frmt', {'svg', 'mat'});

end

function plot_fr_fit_unit(v, t, tPert, grpLbls)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot FR-fit unit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Hard-coded parameters
uIdx = 28;
uIdx = 25;
% uIdx = 24;
uIdx = 18;
uIdx = 20;    
uIdx = 40;

% Combine all frr structures from control group
frr = catfields([v{1}(:).frr], 1, true);
uGood = frr.mdl.uRcv;
mdlName = frr.mdlName{1};

% Set graphics
[hFig, hAx] = plot_set('Firing Rate (Hz)');
clr = mcu_clr;
hold on;

% Plot raw firing rate
fr = frr.fr(uIdx, :);
plot(t, fr, 'Color', [0.5, 0.5, 0.5], 'LineWidth', 0.5,...
    'HandleVisibility', 'off');

% Plot fit curve (only recovery portion) with legend entry
fitCurve = frr.frMdl(uIdx, :);
[~, uRcv] = min(fitCurve);
txtMdl = [upper(mdlName{uIdx}(1)), mdlName{uIdx}(2 : end)];
plot(t(uRcv:end), fitCurve(uRcv:end), '-', ...
    'Color', 'k', 'LineWidth', 2.5, ...
    'DisplayName', sprintf('Model: "%s."', txtMdl));

% Add asterisk at recovery time
uRcv = round(frr.mdl.rcvTime(uIdx) / frr.info.binSize(1)) + frr.mdl.idxTrough(uIdx);
plot(t(uRcv), fitCurve(uRcv), 'o', 'MarkerEdgeColor', 'none',...
    'MarkerFaceColor', 'k', ...
    'MarkerSize', 12, 'LineWidth', 2, ...
    'DisplayName', sprintf('Rcv. Time (hr): %.1f ',...
    round(frr.mdl.rcvTime(uIdx) * 3 / 60 / 60)));

% % Add asterisk at recovery time (model-free)
% uRcv = round(frr.mdlF.rcvTime(uIdx) / frr.info.binSize(1)) + frr.mdlF.idxTrough(uIdx);
% plot(t(uRcv), fr(uRcv), 'diamond', 'MarkerEdgeColor', 'none',...
%     'MarkerFaceColor', 'k', ...
%     'MarkerSize', 12, 'LineWidth', 2, ...
%     'DisplayName', sprintf('Rcv. Time: %d min',...
%     round(frr.mdlF.rcvTime(uIdx) * 3 / 60)));

% Add vertical dashed line for recovery gain
clrLn = bone(6);
tGain = t(end) * 0.9;
frTrough = frr.mdl.frTrough(uIdx);
frSs = frr.mdl.frSs(uIdx);
plot([tGain, tGain], [frTrough, frSs], ':', 'Color', clrLn(4, :), 'LineWidth', 3, ...
    'DisplayName', sprintf('Rcv. Gain: %.2f', frr.mdl.rcvGain(uIdx)));

% Add vertical dashed line for perturbation depth
tTrough = t(frr.mdl.idxTrough(uIdx));
frBsl = frr.mdl.frBsl(uIdx);
plot([tTrough, tTrough], [frTrough, frBsl], '--', 'Color', clrLn(3, :), 'LineWidth', 3, ...
    'DisplayName', sprintf('Pert. Depth: %.2f', frr.mdl.pertDepth(uIdx)));

% Fill area of spk deficit
idxPert = frr.info.idxPert(1);
tRcv = t(idxPert:end);
frRcv = fr(idxPert:end);
bslLine = ones(size(tRcv)) * frr.mdl.frBsl(uIdx);
tRcv = tRcv(:)';
frRcv = frRcv(:)';
bslLine = bslLine(:)';
fill(hAx, [tRcv, fliplr(tRcv)], [bslLine, fliplr(frRcv)],...
    [0.5, 0.5, 0.5], 'FaceAlpha', 0.3, 'EdgeColor', 'none',...
    'DisplayName', sprintf('Spk. Deficit: %.2f', frr.mdl.spkDfct(uIdx)));

% Formatting
plot_axSize('hFig', hFig, 'szOnly', false, 'axShape', 'wide',...
    'axHeight', 300, 'flgPos', true);
hLgd = legend('interpreter', 'tex');
hLgd.Location = "southeast";
hLgd.Location = "northwest";
hLgd.Position(2) = hLgd.Position(2) - 0.03;
hLgd.Position(1) = hAx.Position(1);
hLgd.Box = 'on';

% Mark perturbation
plot_pert(hAx, tPert, clr.bac);

% Save
% fname = ['MEA ~ fitUnit_', num2str(uIdx)];
% lme_save('hFig', hFig, 'fname', fname, 'frmt', {'svg', 'mat'});


end

function plot_pert(hAx, tPert, clrPert)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot perturbation line and label
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yLimit = ylim;
yText = yLimit(2);
yLine = yLimit(2) - 0.05 * (yLimit(2) - yLimit(1));
plot([tPert, 24], [yLine, yLine], 'Color', clrPert,...
    'LineWidth', 4, 'HandleVisibility', 'off');
text(tPert, yText, 'Baclofen', 'FontName', 'Arial', 'FontSize', 16, ...
    'Color', clrPert, 'HorizontalAlignment', 'left');
end

function [hFig, hAx] = plot_set(yLabel)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set axis formatting and create figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[hFig, hAx] = plot_axSize('szOnly', false, 'axShape', 'wide', 'axHeight', 300);
xlabel('Time (Hr)');
ylabel(yLabel);
xlim([-4, 22]);
xticks([0 : 6 : 24]);
end


