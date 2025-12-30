function test_prc_mcmc()
% TEST_PRC_MCMC Verification of MCMC implementation against original.

    % 1. Generate Synthetic Data
    nUnits = 20;
    T = 100; % seconds
    rate = 5; % Hz
    
    fprintf('Generating synthetic spikes (N=%d, T=%ds, Rate=%dHz)...\n', nUnits, T, rate);
    spktimes = cell(nUnits, 1);
    for i = 1:nUnits
        % Create Poisson spike trains
        spktimes{i} = sort(rand(round(T * rate * 1.5), 1) * T);
        spktimes{i} = spktimes{i}(spktimes{i} <= T);
    end
    
    % Shared parameters
    params = {...
        'winLim', [0, T], ...
        'binSize', 0.005, ... % 5ms bins
        'nShuffles', 200, ... % Lower for speed
        'flgSave', false, ...
        'winStpr', 0.5 ...
    };

    % 2. Run Original Implementation
    fprintf('\nRunning ORIGINAL prc_calc...\n');
    tic;
    prc_old = prc_calc_old(spktimes, params{:});
    t_old = toc;
    fprintf('Original done in %.2fs\n', t_old);

    % 3. Run New Implementation
    fprintf('\nRunning NEW prc_calc (MCMC)...\n');
    try
        tic;
        prc_new = prc_calc(spktimes, params{:});
        t_new = toc;
        fprintf('New done in %.2fs\n', t_new);
    catch ME
        fprintf('Error running NEW prc_calc: %s\n', ME.message);
        rethrow(ME);
    end

    % 4. Compare Results
    fprintf('\n--- COMPARISON ---\n');
    
    % A. Real PRC (should be identical, as logic didn't change for data)
    diff_real = max(abs(prc_old.prc0 - prc_new.prc0), [], 'omitnan');
    fprintf('Max diff in Real PRC0: %.6e (Should be 0)\n', diff_real);
    assert(diff_real < 1e-9, 'Real PRC calculation changed!');

    % B. Shuffled Stats
    mu_old = mean(prc_old.prc0_shfl, 2, 'omitnan');
    mu_new = mean(prc_new.prc0_shfl, 2, 'omitnan');
    
    sd_old = std(prc_old.prc0_shfl, [], 2, 'omitnan');
    sd_new = std(prc_new.prc0_shfl, [], 2, 'omitnan');
    
    % Comparison of means (across shuffles)
    err_mu = mean(abs(mu_old - mu_new), 'omitnan');
    fprintf('Mean absolute diff in Shuffled Mean: %.6e\n', err_mu);
    
    % Comparison of stds (across shuffles)
    err_sd = mean(abs(sd_old - sd_new), 'omitnan');
    fprintf('Mean absolute diff in Shuffled Std:  %.6e\n', err_sd);
    
    % Distribution check (KS test for first good unit)
    uTest = prc_old.info.uGood(1);
    [h, p] = kstest2(prc_old.prc0_shfl(uTest, :), prc_new.prc0_shfl(uTest, :));
    fprintf('KS Test for Unit %d: h=%d, p=%.4f (h=0 means distributions are same)\n', uTest, h, p);

    % C. Performance
    fprintf('\nSpeedup: %.2fx\n', t_old / t_new);
    
    if err_mu < 0.5 && err_sd < 0.5
        fprintf('\n[PASS] MCMC implementation gives similar null distribution stats.\n');
    else
        fprintf('\n[WARN] MCMC stats differ significantly.\n');
    end

end
