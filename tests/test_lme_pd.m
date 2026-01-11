% Test LME Partial Dependence (Marginal Effects)
% Verifies lme_pd functionality for 1D and 2D interaction plots.

root = 'd:\Code\slutsky_ECInVivo';
addpath(fullfile(root, 'lme'));

fprintf('---------------------------------------------------\n');
fprintf('TEST: LME Partial Dependence (lme_pd)\n');
fprintf('---------------------------------------------------\n');

%% 1. Create Dummy Data
% Scenario: frSs = frTrough + fr + pBspk*Group + (1|Name)
n = 200;
rng(42);

tbl = table();
tbl.Name = repmat({'M1'; 'M2'; 'M3'; 'M4'; 'M5'}, n/5, 1);
tbl.Group = repmat({'Ctrl'; 'KO'}, n/2, 1);
tbl.fr = randn(n, 1) * 2 + 5;
tbl.frTrough = randn(n, 1) * 1 + 2;
tbl.pBspk = rand(n, 1); % 0 to 1

% Effects
% Group: KO has higher baseline
% pBspk: Positive effect
% Interaction: pBspk effect is stronger in KO
baseline = 10;
eff_fr = 0.5;
eff_pBspk = 5;
eff_Grp = 2; % KO offset
eff_Inter = 5; % Extra slope for KO

isKO = strcmp(tbl.Group, 'KO');
mu = baseline + ...
    eff_fr * tbl.fr + ...
    eff_pBspk * tbl.pBspk + ...
    eff_Grp * isKO + ...
    eff_Inter * (tbl.pBspk .* isKO);

% Random Effect
uNames = unique(tbl.Name);
reMap = containers.Map(uNames, randn(size(uNames)));
reVec = zeros(n, 1);
for i=1:n, reVec(i) = reMap(tbl.Name{i}); end

tbl.frSs = mu + reVec + randn(n, 1)*0.5;

frml = 'frSs ~ frTrough + fr + pBspk * Group + (1|Name)';

%% 2. Fit Model
fprintf('\n--- 2. Fitting Model ---\n');
mdl = lme_fit(tbl, frml, 'dist', 'Normal');
disp(mdl.Coefficients);

%% 3. Test 1D Plot (Continuous)
fprintf('\n--- 3. Testing 1D Plot (pBspk) ---\n');
try
    [pd1, h1] = lme_pd(mdl, {'pBspk'});
    title(h1.CurrentAxes, '1D: pBspk (Marginal)');
    fprintf('[PASS] 1D Plot created.\n');
    disp(head(pd1));
catch ME
    fprintf('[FAIL] 1D Plot error: %s\n', ME.message);
    disp(ME.stack(1));
end

%% 4. Test 2D Interaction (Continuous * Categorical)
fprintf('\n--- 4. Testing 2D Interaction (pBspk * Group) ---\n');
try
    [pd2, h2] = lme_pd(mdl, {'Group', 'pBspk'});
    title(h2.CurrentAxes, '2D: pBspk * Group');
    fprintf('[PASS] 2D Plot created.\n');
    disp(head(pd2));
catch ME
    fprintf('[FAIL] 2D Plot error: %s\n', ME.message);
    disp(ME.stack(1));
end

%% 5. Test 2D Interaction (Continuous * Continuous - pBspk * fr)
% Ensure no crash even if meaningfulness is questionable for this specific fitted model
fprintf('\n--- 5. Testing 2D Interaction (pBspk * fr) ---\n');
try
    % Note: Variable 'fr' must be in vars
    [pd3, h3] = lme_pd(mdl, {'fr', 'pBspk'});
    title(h3.CurrentAxes, '2D: pBspk * fr');
    fprintf('[PASS] 2D Continuous Plot created.\n');
    disp(head(pd3));
catch ME
    fprintf('[FAIL] 2D Continuous Plot error: %s\n', ME.message);
    disp(ME.stack(1));
end

fprintf('\nDone.\n');
