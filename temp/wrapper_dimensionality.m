

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% organize data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% params
thrDim = 0.8;
binsize = 1;

% file
basepath = 'F:\Data\MEA\ketamine\201227_161755';
cd(basepath)
[~, basename] = fileparts(basepath);

% load data
varsFile = ["mea";];
varsName = ["mea";];
v = getSessionVars('basepaths', {basepath}, 'varsFile', varsFile,...
    'varsName', varsName);

% calc firing rate
winStart = 60 * 60;
win = [0, 15 * 60] + winStart;

fr_exp = calc_fr(v.mea.spktimes, 'basepath', basepath,...
    'graphics', false, 'binsize', 60 * 10, 'saveVar', false, 'forceA', true,...
    'smet', 'none', 'winBL', [0, Inf], 'winCalc', [0, Inf]);

fr = calc_fr(v.mea.spktimes, 'basepath', basepath,...
    'graphics', false, 'binsize', binsize, 'saveVar', false, 'forceA', true,...
    'smet', 'none', 'winBL', win, 'winCalc', win);

% firing rate data
fr_mat = fr.strd';
mfr = mean(fr_mat, 2);
nunits = size(fr_mat, 2);

% sort by unit fr
frs = mean(fr_mat, 1);
[~, frIdx] = sort(frs, 'descend');
fr_mat = fr_mat(:, frIdx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calc dimensionality
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% dimensionality
dims = calc_dimensionality(fr_mat, thrDim)

% remove unit fr
fr_mat2 = fr_mat - frs;
fr_mat2(fr_mat2 < 0) = 0;
fr_mat2 = zeros(size(fr_mat));
for iunit = 1 : nunits
    if any(fr_mat(:, iunit))
        fr_mat2(:, iunit) = bz_NormToRange(fr_mat(:, iunit), [0, 1]);
    end
end

% fr_mat2 = [bz_NormToRange(fr_mat, [0, 1])];
dims2 = calc_dimensionality(fr_mat2, thrDim)

% remove the MFR trajectory from each neuron's firing rate
fr_mat3 = zeros(size(fr_mat));
for iunit = 1 : nunits

    % linear regression to find the best scalar
    b = regress(fr_mat2(:, iunit), [mfr, ones(size(mfr))]);

    % subtract the influence of MFR
    fr_mat3(:, iunit) = fr_mat2(:, iunit) - (mfr * b(1));

end
for iunit = 1 : nunits
    if any(fr_mat2(:, iunit))
        fr_mat3(:, iunit) = bz_NormToRange(fr_mat2(:, iunit), [0, 1]);
    end
end


dims3 = calc_dimensionality(fr_mat3, thrDim)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% open figure
setMatlabGraphics(true)
fh = figure;
set(fh, 'WindowState', 'maximized');
tlayout = [3, 3];
th = tiledlayout(tlayout(1), tlayout(2));
th.TileSpacing = 'tight';
th.Padding = 'none';
title(th, basename, 'interpreter', 'none', 'FontSize', 20)
set(fh, 'DefaultAxesFontSize', 16);

% spectrogram
axh1 = nexttile(th, 1, [1, 3]); cla; hold on
plot(fr_exp.tstamps / 60 / 60, fr_exp.strd)
axis tight
yLimit = ylim;
plot(win / 60 / 60, [yLimit(2) yLimit(2)], 'Color', 'k', 'LineWidth', 5)

% map of raw FRs 
axh2 = nexttile(th, tlayout(2) + 1, [1, 1]); cla; hold on
coff = [0, prctile(fr_mat(:), 98)];
PlotColorMap(fr_mat', 1, 'bar','on', 'cutoffs', coff, 'x', fr.tstamps / 60);
ylabel('Unit No.')
title('Raw Firing Rate');
xlabel('Time (m)')

% map normalized FRs
axh3 = nexttile(th, tlayout(2) + 2, [1, 1]); cla; hold on
coff = [0, prctile(fr_mat2(:), 98)];
PlotColorMap(fr_mat2', 1, 'bar','on', 'cutoffs', coff, 'x', fr.tstamps / 60);
ylabel('Unit No.')
title('MFR-Adjusted Firing Rate');
xlabel('Time (m)')

% map of MFR-adjusted FRs
axh4 = nexttile(th, tlayout(2) + 3, [1, 1]); cla; hold on
coff = [0, prctile(fr_mat3(:), 98)];
PlotColorMap(fr_mat3', 1, 'bar','on', 'cutoffs', coff, 'x', fr.tstamps / 60);
ylabel('Unit No.')
title('MFR-Adjusted Firing Rate');
xlabel('Time (m)')

axh5 = nexttile(th, tlayout(2) * 2 + 1, [1, 1]); cla; hold on
plot(fr.tstamps / 60, mean(fr_mat, 2), 'Color', 'k', 'LineWidth', 3)
ylabel('MFR (Hz)')

axh6 = nexttile(th, tlayout(2) * 2 + 2, [1, 1]); cla; hold on
plot(fr.tstamps / 60, mean(fr_mat2, 2), 'Color', 'k', 'LineWidth', 3)
ylabel('MFR (Hz)')

axh7 = nexttile(th, tlayout(2) * 2 + 3, [1, 1]); cla; hold on
plot(fr.tstamps / 60, mean(fr_mat3, 2), 'Color', 'k', 'LineWidth', 3)
ylabel('MFR (Hz)')

linkaxes([axh2, axh5], 'x')
linkaxes([axh3, axh6], 'x')
linkaxes([axh4, axh7], 'x')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sub-functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function dims = calc_dimensionality(mat, thr)

% apply pca to mat
[~, pc, latent, ~, expl] = pca(mat);

% calc explained variance
exp_var = cumsum(latent) / sum(latent);
dims = find(exp_var >= thr, 1, 'first');

fprintf('Number of components needed to explain 80%% of variance: %d\n', dims);

end



































%
%
%
% session = CE_sessionTemplate(pwd, 'viaGUI', false,...
%     'forceDef', true, 'forceL', true, 'saveVar', true);
% nchans = session.extracellular.nChannels;
% fs = session.extracellular.sr;
% spkgrp = session.extracellular.spikeGroups.channels;
% [~, basename] = fileparts(basepath);
%
%
% % load data
% varsFile = ["session"; "datInfo"; "units";];
% varsName = ["session"; "datInfo"; "units"];
% v = getSessionVars('basepaths', {basepath}, 'varsFile', varsFile,...
%     'varsName', varsName);
%
% % firing rate
% load([basename, '.spikes.cellinfo.mat'])
%
% win = [0, 60 * 15];
% fr = calc_fr(spikes.times, 'basepath', basepath,...
%     'graphics', false, 'binsize', 0.1, 'saveVar', false, 'forceA', true,...
%     'smet', 'none', 'winBL', win, 'winCalc', win);
%
% fh = figure;
% plot(fr.tstamps / 60, fr.strd)
% % plot fr vs. time
%
%
% % Example data: fr_matrix is a matrix where rows are time bins and columns are neurons
% % Random data for illustration
% fr_matrix = fr.strd(units.clean(1, :), :)'; % 100 time bins, 150 neurons
% % fr_matrix = fr.strd';
%
% % Step 4: Apply PCA to the adjusted firing rates
% [~, pc, latent, ~, expl] = pca(fr_matrix);
%
% % Step 5: Determine the number of components needed to explain 80% of the variance
% explained_variance = cumsum(latent) / sum(latent);
% explained_variance = cumsum(expl) / 100;
% num_components = find(explained_variance >= 0.8, 1, 'first');
%
% fprintf('Number of components needed to explain 80%% of variance: %d\n', num_components);
%
% % Step 1: Compute the Mean Firing Rate (MFR) for each time bin
% mfr = mean(fr_matrix, 2); % Mean across columns (neurons)
%
% % Step 2: Initialize matrix to store adjusted firing rates
% adjusted_fr_matrix = zeros(size(fr_matrix));
%
% % Step 3: Remove the influence of MFR from each neuron's firing rate
% for i = 1:size(fr_matrix, 2) % Iterate through each neuron
%     % Linear regression to find the best scalar
%     b = regress(fr_matrix(:, i), [mfr, ones(size(mfr))]); % Including intercept
%     % Subtract the influence of MFR
%     adjusted_fr_matrix(:, i) = fr_matrix(:, i) - (mfr * b(1));
% end
%
% % Step 4: Apply PCA to the adjusted firing rates
% [coeff, score, latent] = pca(adjusted_fr_matrix);
%
% % Step 5: Determine the number of components needed to explain 80% of the variance
% explained_variance = cumsum(latent) / sum(latent);
% num_components = find(explained_variance >= 0.8, 1, 'first');
%
% % Display results
% fprintf('Number of components needed to explain 80%% of variance: %d\n', num_components);