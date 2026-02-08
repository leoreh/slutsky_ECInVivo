%% Debug QQ Plot Cutoff
% This script visualizes the difference between filtering by empirical
% percentile vs theoretical log-normal distribution properties.

rng(42);
% Generate log-normal data: mu=1, sigma=1 (log scale)
% 1000 points
fr = lognrnd(1, 1, 1000, 1);

% Introduce a "floor" or minimum value often seen in real data
% e.g. min firing rate detectable is 0.01 Hz
fr(fr < 0.01) = 0.01; 

cutoff_z = -4;
cutoff_p = normcdf(cutoff_z); % ~3.16e-5

% 1. Empirical Percentile (Current Logic)
cutoff_Hz_emp = prctile(fr, cutoff_p * 100);

% 2. Theoretical Log-Normal Cutoff
log_fr = log(fr);
mu_hat = mean(log_fr);
sigma_hat = std(log_fr);

% Theoretical log-value at Z = -4
log_cutoff_theo = mu_hat + (cutoff_z * sigma_hat);
cutoff_Hz_theo = exp(log_cutoff_theo);

fprintf('Cutoff Z: %.2f (p = %.2e)\n', cutoff_z, cutoff_p);
fprintf('Empirical Cutoff (prctile): %.5f Hz\n', cutoff_Hz_emp);
fprintf('Theoretical Cutoff (mu+Zs): %.5f Hz\n', cutoff_Hz_theo);

% Count removed
n_total = length(fr);
n_removed_emp = sum(fr < cutoff_Hz_emp);
n_removed_theo = sum(fr < cutoff_Hz_theo);

fprintf('Units removed (Empirical):  %d / %d\n', n_removed_emp, n_total);
fprintf('Units removed (Theoretical): %d / %d\n', n_removed_theo, n_total);

% Plot
figure;
h = qqplot(log(fr));
hold on;
yline(log(cutoff_Hz_emp), 'r-', 'Empirical (prctile)');
yline(log_cutoff_theo, 'b--', 'Theoretical (mu+Zs)');
legend('Data', 'Fit', 'Empirical Cutoff', 'Theoretical Cutoff', 'Location', 'best');
title('Debug QQ Cutoff');
grid on;
