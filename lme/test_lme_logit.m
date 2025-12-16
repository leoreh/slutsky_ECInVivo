% Test script for lme_compareDists Logit-Normal addition

% Create dummy data
rng(42);
n = 100;
tbl = table();
tbl.Group = randi(2, n, 1); % Random factor

% 1. Beta distributed data (ideal for Logit-Normal)
% Beta(2,2) is bell-ish on [0,1]
tbl.yBeta = betarnd(2, 2, n, 1);

% 2. Gamma distributed data (positive, skewed, >1 possible)
tbl.yGamma = gamrnd(2, 2, n, 1);

% 3. Data with zeros (for validation check)
tbl.yZeros = tbl.yBeta;
tbl.yZeros(1:5) = 0;

disp('--- Testing Beta Data (0,1) ---');
statsBeta = lme_compareDists(tbl, 'yBeta ~ 1 + (1|Group)');
disp(statsBeta);

disp('--- Testing Gamma Data (should have Logit-Normal but maybe poor fit if > 1 clipped) ---');
try
    statsGamma = lme_compareDists(tbl, 'yGamma ~ 1 + (1|Group)');
    disp(statsGamma);
catch ME
    disp(['Gamma Error: ' ME.message]);
end

disp('--- Testing Zero Data (Validation) ---');
try
    statsZeros = lme_compareDists(tbl, 'yZeros ~ 1 + (1|Group)');
    disp(statsZeros);
catch ME
    disp(['Zero Error: ' ME.message]);
end
