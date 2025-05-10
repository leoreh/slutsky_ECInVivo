function [boutLbls, boutStats] = bouts_separate(bouts, varargin)
%BOUTS_SEPARATE Classifies bout durations into spike and slab categories using different methods.
%
% Overview:
% This function takes a matrix of bout start and end times and classifies each bout
% into one of two categories: "spike" (shorter durations) or "slab" (longer
% durations). Two primary methods for this classification are supported, selectable
% via the 'sepMet' parameter:
%   1. 'gmm': Gaussian Mixture Model fitting on log-transformed durations.
%   2. 'mcshane': Threshold-based separation with Gamma fitting for slab durations,
%      adapted from McShane et al. (2010, Journal of Neuroscience Methods).
%
% Method Details:
%
%   'gmm' Method:
%   - Operates on log-transformed bout durations to better handle skewed distributions.
%   - Fits a 2-component Gaussian Mixture Model (GMM) to these log-durations.
%   - The two components are assumed to represent "spike" (smaller mean in log-domain)
%     and "slab" (larger mean in log-domain) populations.
%   - Bouts are classified based on the posterior probability of belonging to each
%     component; a bout is assigned to the component with the higher posterior probability.
%   - `boutDist` output for GMM:
%       .props: [1x2] Sorted mixing proportions of the components [spike_prop; slab_prop].
%       .fit: (gmdistribution object) The fitted GMM object.
%       .thrSlab: NaN (not applicable to GMM method).
%
%   'mcshane' Method:
%   - Inspired by McShane et al. (2010), this method uses a duration threshold.
%   - Bouts with duration <= `thrSlab` (default 40 units, applied to original durations)
%     are classified as "spike".
%   - Bouts with duration > `thrSlab` are classified as "slab".
%   - A Gamma distribution is fitted to the durations of the "slab" bouts.
%   - `boutDist` output for McShane:
%       .props: [1x2] Proportions of spike and slab bouts based on the threshold [spike_prop; slab_prop].
%       .fit: [1x3 vector] Parameters of the Gamma distribution fitted to slab durations
%             [shape_a, scale_b, mean_gamma]. Set to [NaN NaN NaN] if fit fails or no slab bouts.
%       .thrSlab: (scalar) The threshold used for classification.
%
% Common `boutDist` Output Fields (for both methods, based on original durations):
%   .avg: [1x2] Mean duration of [spike_bouts, slab_bouts].
%   .var: [1x2] Variance of duration of [spike_bouts, slab_bouts].
%   .durTotal: [1x2] Total duration of [spike_bouts, slab_bouts].
%
% General Considerations:
%   - Individual Variability: For analyzing multiple subjects, it is generally recommended
%     to call this function separately for each subject\'s data to account for
%     inter-individual differences in bout duration distributions.
%   - Input Units: Ensure consistency in the time units of the input `bouts` matrix,
%     as this will affect duration calculations and the interpretation of `thrSlab`.
%
% INPUT:
%   bouts           (numeric matrix) Nx2 matrix where N is the number of
%                   bouts. Each row is [start_time, end_time] of a
%                   bout in seconds.
%   varargin        Optional name-value pairs:
%     'flgGraphics' (logical) If true, a figure visualizing the classification
%                     is generated. Defaults to false.
%     'sepMet'      (char/string) Separation method: 'gmm' (default) or 'mcshane'.
%
% OUTPUT:
%   boutType        (logical vector) Nx1 logical vector.
%                   True (1) if a bout is classified as "slab".
%                   False (0) if a bout is classified as "spike".
%   boutDist        (struct) Structure containing parameters of the
%                   fitted populations, depending on the method of
%                   separation (see above for details).
%
% Dependencies:
%   Statistics and Machine Learning Toolbox (for `fitgmdist`, `posterior`,
%   `normpdf`, `gampdf`, `gamfit`)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Parser
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addRequired(p, 'bouts', @(x) isnumeric(x) && (size(x,2) == 2 || isempty(x)));
addOptional(p, 'flgGraphics', false, @islogical);
addOptional(p, 'sepMet', 'GMM', @(x) (ischar(x) || isstring(x)) && ismember(lower(x), {'gmm', 'mcshane'}));

parse(p, bouts, varargin{:});

flgGraphics             = p.Results.flgGraphics;
sepMet                  = lower(p.Results.sepMet); % Use lowercase for reliable comparison

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate bout durations
durs = bouts(:,2) - bouts(:,1);
dursLog = log(durs);

nBouts = size(bouts, 1);
boutLbls = false(nBouts, 1); % Default to false (spike)

% Initialize boutDist structure
boutStats = struct();
boutStats.avg = [NaN NaN];
boutStats.var = [NaN NaN];
boutStats.nBouts = [NaN NaN];
boutStats.dur = [NaN NaN];
boutStats.durTotal = [NaN NaN];
boutStats.props = [NaN NaN]; 
boutStats.fit = [];        
boutStats.thrSlab = NaN;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Classification
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch sepMet
    case 'gmm'
        
        % Fit a 2-component Gaussian Mixture Model
        opts = statset('MaxIter', 1000, 'TolFun', 1e-6);
        gm = fitgmdist(dursLog(:), 2, 'Options', opts,...
            'CovarianceType', 'diagonal', 'SharedCovariance', false);

        % Sort parameters by mean (spike first, then slab)
        [logMu, idxSort] = sort(squeeze(gm.mu));
        
        logSigma = squeeze(gm.Sigma);
        logSigma = logSigma(idxSort); % [spike_var_log; slab_var_log]
        
        current_props = gm.ComponentProportion';
        boutStats.props = current_props(idxSort)'; % Ensure 1x2 [spike_prop, slab_prop]
        
        boutStats.fit = gm;

        post_probs = posterior(gm, dursLog(:));
        [~, idxType] = max(post_probs, [], 2);

        % boutType is true (1) if assigned to slab component.
        % idxSort(2) is the original index of the component identified as slab (larger mean).
        boutLbls = (idxType == idxSort(2));

    case 'mcshane'
        thrSlab = 40; % Duration in original units (e.g., seconds or epochs)
        boutStats.thrSlab = thrSlab;

        % Classify based on thrSlab
        spkIdx = (durs <= thrSlab);
        boutLbls = ~spkIdx; % True for slab, False for spike

        nSpikes = sum(spkIdx);
        nSlabs = nBouts - nSpikes;
        
        boutStats.props = [nSpikes/nBouts, nSlabs/nBouts]; % [spike_prop, slab_prop]

        % Fit gamma distribution to slab durations (if any)
        slab_durs_mc = durs(boutLbls);
        if ~isempty(slab_durs_mc) && numel(slab_durs_mc) >= 2 % gamfit needs at least 2 points
            params_gam = gamfit(slab_durs_mc);
            boutStats.fit = [params_gam, params_gam(1) * params_gam(2)]; % [shape, scale, mean_gamma]
        else
            boutStats.fit = [NaN, NaN, NaN]; % Indicate fit failure or no slab data
        end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate Common Descriptive Statistics (Original Durations)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

spkDurs = durs(~boutLbls);
slabDurs = durs(boutLbls);

boutStats.avg = [mean(spkDurs), mean(slabDurs)]; 
boutStats.var = [var(spkDurs), var(slabDurs)];   
boutStats.dur = cell2padmat({spkDurs, slabDurs}, 2);
boutStats.durTotal = [sum(spkDurs), sum(slabDurs)]; 
boutStats.nBouts = [length(spkDurs), length(slabDurs)]; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Visualization 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flgGraphics
    figure; % Create a new figure for the plot
    
    if strcmp(sepMet, 'gmm')
        hold on;
        histogram(dursLog, 'Normalization', 'pdf', 'DisplayName', 'Log Durations');
        
        x_grid = linspace(min(dursLog), max(dursLog), 200)';
        
        % PDF of the spike component (component 1 in sorted boutDist)
        if ~isempty(logMu) && ~isempty(logSigma) && ~any(isnan(boutStats.props)) % Check if boutDist GMM fields are populated
            pdf_spike_comp = normpdf(x_grid, logMu(1), sqrt(logSigma(1)));
            plot(x_grid, boutStats.props(1) * pdf_spike_comp, 'r-', 'LineWidth', 2,...
                 'DisplayName', sprintf('Spike (Prop: %.2f)', boutStats.props(1)));
        
            % PDF of the slab component (component 2 in sorted boutDist)
            pdf_slab_comp = normpdf(x_grid, logMu(2), sqrt(logSigma(2)));
            plot(x_grid, boutStats.props(2) * pdf_slab_comp, 'g-', 'LineWidth', 2,...
                 'DisplayName', sprintf('Slab (Prop: %.2f)', boutStats.props(2)));
        
            % PDF of the overall mixture
            if isa(boutStats.fit, 'gmdistribution')
                pdf_mixture = pdf(boutStats.fit, x_grid);
                plot(x_grid, pdf_mixture, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Overall Mixture');
            end
        else
            disp('GMM parameters not available for plotting.');
        end
        
        title('GMM Classification of Bout Durations (Log-Domain)');
        xlabel('Log(Duration)');
        ylabel('Probability Density');
        legend('show', 'Location', 'best');
        grid on;
        hold off;

    elseif strcmp(sepMet, 'mcshane')
        hold on;
        maxEpochVal = max(durs);
        edges = logspace(0, log10(maxEpochVal + 0.5), 200);
        histogram(durs, edges, 'Normalization', 'probability',...
            'DisplayName', 'All Bouts', 'FaceAlpha', 1, 'EdgeColor', 'none');

        x_grid_slab = linspace(boutStats.thrSlab + 0.5, maxEpochVal + 0.5, 200);
        pdf_slab_fitted = gampdf(x_grid_slab, boutStats.fit(1), boutStats.fit(2));
        plot(x_grid_slab, boutStats.props(2) * pdf_slab_fitted, 'g-', 'LineWidth', 2, ...
            'DisplayName', sprintf('Slab Gamma Fit (Prop: %.2f)', boutStats.props(2)));

        xlim_upper = maxEpochVal + 0.5;
        xlim([0.5, xlim_upper]);
        title(sprintf('McShane Classification (Threshold = %.1f)', boutStats.thrSlab));
        xlabel('Duration');
        ylabel('Probability');
        legend('show', 'Location', 'best');
        grid on;
        hold off;
    end
end

end
