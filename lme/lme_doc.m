
%% ========================================================================
%  NOTE: EFFECTS
%  ========================================================================
% Linear Mixed-Effects (LME) models are used for hierarchically structured
% data (e.g., repeated measures within subjects) to account for
% non-independence of observations. They partition variance into fixed and
% random effects. This function fits models using Restricted Maximum
% Likelihood (REML) by default, as opposed to Maximum Likelihood (ML),
% which provides unbiased estimates of variance components.
%
% Fixed Effects', denoted 'beta', represent the average effects of
% predictor variables across the entire population being studied. These are
% usually the primary focus of the analysis.
%
% 'Random Effects', denoted 'b', capture the variability *between* the
% groups or subjects. They model how individual subjects systematically
% deviate from the average population. For instance, a 'random intercept'
% allows each subject to have their own unique baseline response level
% (e.g., specified as `(1 | Subject)` in the formula). A 'random slope'
% allows the effect of a predictor, like time, to differ across subjects
% (e.g., `(Time | Subject)`). By incorporating random effects, the model
% properly accounts for the data's dependency structure when evaluating the
% fixed effects.
%
% 'Residual Error', denoted as 'epsilon', accounts for the
% remaining, unexplained variability or random noise in the data that isn't
% captured by either the fixed or the random effects structure.
%  ========================================================================


%% ========================================================================
%  NOTE: G/LME
%  ========================================================================
% The choice between using LME and GLME is determined exclusively by
% the dist of the response variable. LME is used when the
% response variable is continuous and follows a normal dist, or can
% be transformed to approximate one. GLME (Generalized Linear
% Mixed-Effects Model) is required when the response variable does not
% follow a normal dist.

% The decision to transform predictor variables is a separate consideration
% from the choice of model and is applicable to both LME and GLME.
% Predictors should be log-transformed when their dist is highly
% skewed. The purpose of this transformation is twofold: it makes the
% predictor's dist more symmetric, and more importantly, it can
% linearize a non-linear relationship between the predictor and the
% response. By converting an exponential or power-law relationship into a
% linear one, the log-transformed predictor better satisfies the
% fundamental assumption of linear models. This practice reduces the
% disproportionate influence of extreme data points and often leads to a
% more accurate and stable model fit.
%
% Z-scoring predictors is performed not to address dist shape but
% to standardize variables onto a common scale. This is particularly useful
% when predictors are measured in different units (e.g., firing rate in Hz,
% time in seconds). After z-scoring, the predictors model coefficients are
% directly comparable, as a larger coefficient now unambiguously indicates
% a stronger effect per standard deviation, regardless of the original
% units.

% TLDR: Z-score and log-transform contineous variables only when they are
% predictors. As response variables, adjust the model to their
% dist.
%  ========================================================================


%% ========================================================================
%  NOTE: DUMMY CODING & COEFFICIENTS
%  ========================================================================
% MATLAB's `fitlme` and `fitglme` functions use 'Dummy Coding' by default
% for categorical predictors. (also known as 'treatment contrasts'). This
% means one level of each categorical factor is chosen as the 'reference
% level' (often the first alphabetically or numerically). All other levels
% of that factor are then compared directly to this reference level.

% The '(Intercept)' coefficient represents the estimated average value of
% the response variable under a specific baseline condition: when all
% categorical predictors in the model are at their reference levels, and
% simultaneously, all continuous predictors included in the model formula
% are equal to zero. A statistically significant intercept test suggests
% this baseline mean is different from zero.
%
% A 'Main Effect' coefficient, such as one named 'FactorA_LevelX',
% estimates the *difference* in the average response between 'LevelX' of
% 'FactorA' and the reference level of 'FactorA'. However, a critical point
% arises when interactions are present in the model: this comparison is
% valid *only* when all other factors that interact with 'FactorA' are held
% at their respective reference levels. Therefore, in the presence of
% interactions, this coefficient does not represent the overall main effect
% averaged across the levels of other factors.
%
% An 'Interaction' Coefficient, for example
% 'FactorA_LevelX:FactorB_LevelY', assesses how the effect of one factor
% changes depending on the level of another. It quantifies the *additional
% difference* in the response associated with FactorA being at LevelX (vs.
% its reference) when FactorB is specifically at LevelY (compared to when
% FactorB is at its reference level). It essentially captures a "difference
% of differences". A significant interaction term indicates that the effect
% of FactorA is not constant but depends on the specific level of FactorB
% (and vice-versa).
%  ========================================================================


%% ========================================================================
%  NOTE: CONTRASTS
%  ========================================================================
% This function facilitates hypothesis testing about specific comparisons
% or combinations of the fixed-effect coefficients (the betas).
% Mathematically, each test evaluates a hypothesis like H0: H * beta = 0,
% where 'H' is a row vector that defines the linear combination of
% coefficients forming the contrast of interest. The results table
% (`lme_results`) categorizes these tests into distinct types.
%
% The "Coefficient" type refers to the standard tests for individual model
% coefficients equalling zero. The results for these are taken directly
% from the standard output of the LME/GLME model. No separate contrast
% vector 'H' is constructed or tested using `coefTest` here and these
% standard coefficient tests are *not* included in the multiple comparison
% correction procedures; the `pAdj` column will always be NaN for
% rows of Type "Coefficient".
%
% The "Simple Effect" type is generated specifically when two factors
% interact. It tests the effect of one factor at a single, fixed level of
% the other factor. For instance, it might test the difference between
% LevelX and the reference level of FactorA, specifically *only* when
% FactorB is at LevelY. These tests require constructing a specific
% contrast vector 'H' that typically combines the main effect coefficient
% for FactorA_LevelX with the relevant interaction coefficient (e.g.,
% FactorA_LevelX:FactorB_LevelY). The function runs `coefTest` using this H
% vector. The Estimate (calculated as H * beta_hat) and its Standard Error
% (SE, calculated as sqrt(H * Cov(beta_hat) * H')) are computed. These
% simple effect tests undergo multiple comparison correction if they are
% selected for inclusion in the final results table (via the 'contrasts'
% option). The correction is applied collectively across all selected
% simple and marginal effects.
%
% The "Marginal Effect" type is also generated for two-way interactions. It
% tests the effect of one factor averaged across all the levels of the
% other factor it interacts with. This often aligns more closely with the
% intuitive notion of a "main effect" in the presence of an interaction.
% For example, it might test the difference between LevelX and the
% reference level of FactorA, averaged across all levels of FactorB.
% Constructing the H vector for this involves combining the main effect
% coefficient with a weighted sum of the relevant interaction coefficients,
% where the weights ensure averaging (typically 1 divided by the number of
% levels of the factor being averaged over). Again, `coefTest` is used with
% this H vector, the Estimate and SE are computed, and the H vector is
% stored. Like simple effects, these marginal effect tests undergo multiple
% comparison correction if selected for the output.
%
% All results produced by this function – whether standard coefficients,
% simple effects, or marginal effects – pertain strictly to the *fixed
% effects* part of the LME model. They allow inferences about
% *population averages*, representing overall trends or differences in the
% population from which the sample was drawn. These results do not describe
% effects specific to individual subjects or clusters (often called
% conditional effects). Analyzing subject-specific deviations would require
% examining the estimated random effects ('b') themselves, e.g. using
% `randomEffects(lme)`.
%  ========================================================================


%% ========================================================================
%  NOTE: ANOVA
%  ========================================================================
% The `anova(lme)` method provides F-tests for each term. For a model
% fitted with dummy coding:
%
% Interaction F-test (e.g., for Group:Day): This is a robust test for the
% overall significance of the interaction. It should be examined first.
%
% Main Effect F-tests (e.g., for Group, for Day): These test hypotheses
% about the simple effects at reference levels (due to dummy coding). If
% the interaction is significant, these main effect F-tests are less
% directly interpretable as "averaged effects" and the focus should shift
% to simple/marginal effects that dissect the interaction. If the
% interaction is *not* significant, these F-tests become more interpretable
% as overall main effects.
%
% Note: To obtain Type III F-tests where main effects are tested more akin
% to being "averaged over" other factors in the presence of interactions,
% `fitlme`/`fitglme` would need to be called with `'DummyVarCoding',
% 'effects'`. This function currently uses default dummy coding to simplify
% H-vector construction for simple/marginal effects.
%  ========================================================================


%**************************************************************************
% Estimate +/ SE
%**************************************************************************
% The `Estimate` column provides the model's calculated magnitude for each
% specific term or contrast. For the Intercept, the `Estimate` represents
% the predicted mean value of the response variable when all categorical
% fixed factors are at their designated reference levels and any continuous
% fixed predictors are held at zero. For other coefficients (main effects
% or interactions under dummy coding), the `Estimate` represents the
% calculated difference compared to the relevant reference level(s). For
% instance, a main effect coefficient's Estimate is the difference between
% that factor level and its reference level, specifically evaluated when
% other interacting factors are at their reference levels. An interaction
% coefficient's Estimate quantifies how the effect of one factor changes
% across levels of another. For derived contrasts (like simple or marginal
% effects), the Estimate is the calculated value of the specific linear
% combination of coefficients being tested (e.g., the estimated mean
% difference between two groups at a specific condition). The `SE`
% (Standard Error) accompanies each Estimate and quantifies the precision
% of that estimate; it reflects the expected standard deviation of the
% estimate if the study were repeated many times. The `SE` is crucial for
% assessing statistical significance, as it forms the denominator in the
% t-statistic calculation (t = Estimate / SE) and is used to construct
% confidence intervals around the `Estimate`.
%  ========================================================================


%% ========================================================================
%  NOTE: STATISTIC (T VS F)
%  ========================================================================
% The `Statistic` column contains the test statistic associated with the
% hypothesis test for that row. For rows of type 'Coeff', 'Simple', or
% 'Marginal', this is the t-statistic (Estimate / SE), testing whether the
% specific coefficient or contrast is significantly different from zero.
% For rows of type 'ANOVA', this is the F-statistic obtained from the
% `anova(lme)` function, testing the overall significance of the
% interaction term in the model.
%  ========================================================================


%% ========================================================================
%  NOTE: DEGREES OF FREEDOM (DF)
%  ========================================================================
%  The `DF` column reports the degrees of freedom associated with the test
%  statistic, which reflects the amount of independent information available
%  to estimate the variability.
%
%  Interpretation of the Output:
%     - For t-tests ('Coeff', 'Simple', 'Marginal'): This is a single value
%       representing the denominator degrees of freedom (DF_den).
%     - For F-tests ('ANOVA'): This displays '[DF1, DF2]'. DF1 is the
%       numerator DF (parameters in the effect), and DF2 is the denominator
%       DF (error/residual).
%
%  The Importance of Satterthwaite (for LME):
%     In hierarchical data, observations within the same group (e.g.,
%     subject) are not independent. The default 'Residual' method
%     calculates DF as simply (Total Observations - Number of Parameters).
%     This ignores the correlation structure, drastically overestimating
%     the DF and leading to "optimistic" (false positive) p-values.
%
%     Satterthwaite's approximation corrects for this by estimating the
%     "effective" degrees of freedom based on the variance components. It
%     is essential for unbalanced designs or small sample sizes (e.g., few
%     subjects), providing a more conservative and reliable inference.
%
%  Limitation in GLME (Gamma/Poisson/Log-Normal):
%     MATLAB's `fitglme` does NOT support Satterthwaite approximations
%     because exact finite-sample distributions do not exist for
%     generalized models. For these models, the DF is calculated using the
%     'Residual' method or treated as infinite (Z-tests). Consequently,
%     p-values in GLMEs are based on asymptotic assumptions and should be
%     interpreted with this limitation in mind.
%  ========================================================================


%% ========================================================================
%  NOTE: MULTIPLE COMPARISONS
%  ========================================================================
% Applied only to the family of selected "Simple Effect" and "Marginal
% Effect". Not applied to ANOVA F-tests or individual model coefficients.
%
% To counteract this inflation of the error rate across the set of tests,
% correction methods are applied to the calculated p-values. The
% 'correction' parameter allows choosing between different strategies:
%
% Bonferroni ('bonferroni'): This is the simplest method, controlling the
% Family-Wise Error Rate (FWER) – the probability of making one or more
% Type I errors across all tests performed. It achieves this by multiplying
% each individual p-value by the total number of tests conducted (m) or,
% equivalently, by dividing the target significance level (e.g., 0.05) by
% m. Bonferroni is often overly conservative, significantly reducing
% statistical power and increasing the risk of failing to detect true
% effects (Type II errors).
%
% Holm-Bonferroni ('holm', default): This method also controls the FWER but
% is uniformly more powerful (less conservative) than the standard
% Bonferroni procedure. It works sequentially: p-values are ordered from
% smallest to largest, and the significance threshold is adjusted
% step-by-step. The smallest p-value is compared against alpha/m, the next
% smallest against alpha/(m-1), and so on, stopping when a p-value fails to
% meet its adjusted threshold. This provides the same strong FWER control
% as Bonferroni but with a better chance of detecting true effects.
%
% False Discovery Rate ('fdr', specifically Benjamini-Hochberg): Instead of
% controlling the probability of making any Type I error (FWER), the FDR
% approach controls the expected proportion of rejected null hypotheses
% that are actually false positives. For example, setting FDR control at 5%
% aims to ensure that, among all the effects declared significant, no more
% than 5% are expected to be false discoveries. This method is less
% stringent than FWER control, particularly when a large number of tests
% are performed. Consequently, it offers considerably more power to detect
% true effects, making it suitable for more exploratory analyses where
% controlling the proportion of false findings is deemed acceptable, rather
% than strictly preventing any single false positive.
%  ========================================================================


%% ========================================================================
%  NOTE: FITMETHOD
%  ========================================================================
%  Specifies the objective function maximized to estimate model parameters.
%  The choice depends on whether the goal is Model Selection (comparing
%  different fixed effects) or Parameter Estimation (accurate variances).
%
%  Maximum Likelihood (ML) & Integral Approximation:
%     - LME: 'ML'
%     - GLME: 'Laplace'
%
%     These methods approximate the true marginal likelihood of the data.
%     MANDATORY for comparing models with *different fixed effects* or
%     different distributions using Information Criteria (AIC/BIC).
%     However, variance components (random effects) are biased downwards in
%     finite samples.
%
%  Restricted / Pseudo-Likelihood:
%     - LME: 'REML'  (Restricted Maximum Likelihood)
%     - GLME: 'REMPL' (Restricted Pseudo-Likelihood) or 'MPL'
%
%     These methods estimate parameters based on the likelihood of the
%     residuals (removing fixed effects first) or via linearization.
%     Preferred for final reporting as they provide unbiased estimates of
%     variance components. The resulting likelihoods are NOT comparable
%     across models with different fixed effects. Do NOT use AIC/BIC from
%     these methods.
%  ========================================================================


%% ========================================================================
%  NOTE: HANDLING ZEROS
%  ========================================================================
%  Addresses the constraint that Log-link, Gamma, and Log-Normal models are
%  undefined for y = 0. A common heuristic is transforming y -> log(y + c).
%
%  The "Half-Minimum" Heuristic:
%     Commonly, c is set to min(y > 0) / 2. This creates a floor for the
%     data, allowing log-transformation (or Gamma fitting) while preserving
%     the rank order of observations.
%
%  2. Statistical Implications & Risks:
%     - Sensitivity: The estimated coefficients (elasticities) for small y
%       are highly sensitive to the arbitrary choice of c.
%     - Bias: It introduces a "hard floor" (artificial curvature) in the
%       Mean-Variance relationship, potentially violating model assumptions.
%     - Scale Invariance: The transformation is not scale-invariant (e.g.,
%       adding 0.5 yields different results if y is measured in ms vs sec).
%
%  3. Best Practices:
%     - Sensitivity Analysis: If using a shift, verify results are robust
%       to changes in c (e.g., test c = min(y) and c = min(y)/10).
%     - Modern Alternatives: If zeros are prevalent, prefer distributions
%       that handle zeros naturally (e.g., Tweedie, Poisson/PPML) to avoid
%       arbitrary shifts entirely.
%  ========================================================================

%% ========================================================================
%  NOTE: POOLED VS. GROUP-WISE Z-SCORING
%  ========================================================================
% When standardizing continuous predictors (Z-scoring), the choice between
% pooling all data or Z-scoring within groups (e.g., WT vs. KO)
% fundamentally changes the hypothesis being tested.
%
% 1. Pooled Z-scoring (Standard Approach):
%    Predictors are standardized using the mean and standard deviation of 
%    the entire dataset. This preserves the "Between-Group" variance. If 
%    KO units have a lower baseline firing rate than WT units, pooled 
%    Z-scoring maintains this offset. In a model like `Rcv ~ Group * FR`, 
%    the `Group` coefficient will accurately reflect the difference in 
%    recovery probability between the two genotypes at their natural scales.
%
% 2. Group-wise Z-scoring (Centering Within Context):
%    Predictors are standardized separately for each group. This process 
%    force-centers every group at a mean of 0. Technically, this removes 
%    the main effect of the group from that specific predictor. This is 
%    useful ONLY if you want to study "Relative" effects—for instance, 
%    determining if a unit that is "high-firing for its specific culture" 
%    behaves differently, regardless of the culture's absolute firing rate.
%
% 3. The "Classification" Guardrail:
%    If the `Group` is your response variable (e.g., using FR to predict 
%    if a unit is WT or KO), you must NEVER Z-score within groups. This 
%    constitutes "Data Leakage," as you would be using the category labels 
%    to transform the predictors, effectively hiding the differences the 
%    model is trying to find.
% ========================================================================


%% ========================================================================
%  NOTE: LINK FUNCTIONS
%  ========================================================================
%  A GLM specifies the relationship between the Mean Response (mu) and the
%  linear predictor (X*beta) via a Link Function g(mu) = X*beta.
%
%  While statistical software defaults to "canonical" links for
%  computational simplicity—such as the Reciprocal link for Gamma
%  distributions - these defaults often contradict the physical reality of
%  biological data.

%  The Canonical links for Gamma/InvGauss imply that as X increases, Y
%  decreases hyperbolically (1/X). This is rarely the true biological
%  mechanism. Most positive data (like RT) follows a multiplicative process
%  (e.g., "a 10% increase"), which implies an Exponential trend. This
%  corresponds to the LOG LINK: log(mu) = X*beta.

%  Therefore, the Log link is the scientifically appropriate choice for
%  Gamma, Inverse Gaussian, and Poisson distributions as it models these
%  exponential trends. In contrast, the Normal distribution uses the
%  Identity link by default, representing the standard assumption of
%  additive, linear effects. By enforcing the Log link across the
%  heavy-tailed distributions while keeping Identity for the Normal
%  distribution, we ensure that the AIC comparison properly tests the
%  variance structure (distribution) rather than inadvertently penalizing a
%  plausible distribution for using an implausible hyperbolic trend.
%  ========================================================================


%% ========================================================================
%  NOTE: MODULATION INDEX (MI) VS. LOG-RATIO
%  ========================================================================
%  When quantifying firing rate changes (e.g., Ripple vs. Random periods),
%  the choice of metric significantly impacts the statistical validity
%  of the subsequent Linear Mixed-Effects (LME) analysis.
%
%  1. The Modulation Index (MI):
%     Defined as (FR_a - FR_b) / (FR_a + FR_b). This metric is a standard
%     "normalized difference" in neuroscience because it is intuitive and
%     controls for the absolute firing rate of a neuron. However, it
%     presents challenges for LME models. Because it is strictly bounded
%     between [-1, 1], it violates the assumption of normality, especially
%     near the boundaries. Analyzing MI typically requires rescaling the
%     data to (0, 1) and employing a Logit-Normal model.
%
%  2. The Log-Ratio (LR):
%     Defined as ln(FR_a / FR_b). This is often the statistically preferred
%     alternative to the MI. Unlike the MI, the Log-Ratio is not bounded
%     and is naturally symmetric around zero, which often results in
%     residuals that better approximate a Normal distribution. This
%     symmetry aligns with the "Log Link" philosophy, as it models
%     multiplicative processes on a linear scale.
%
%  3. The GLMM Approach (Raw Rates):
%     The most rigorous method is to avoid index-based metrics entirely. By
%     fitting a Generalized Linear Mixed Model (GLMM) directly to the raw
%     firing rates (e.g., using Gamma or Inverse Gaussian distributions),
%     one preserves the inherent mean-variance relationship of the data
%    . This approach accounts for the "multiplicative" nature
%     of neural spiking while avoiding the information loss and potential
%     biases (like the "retransformation problem") associated with
%     collapsing observations into a single ratio.
%
%  DECISION RULE: Use the Modulation Index (MI) primarily for visualization
%  and consistency with literature. For LME modeling, favor the
%  **Log-Ratio** for simplicity of insference, or a **GLMM on raw rates**
%  (Gamma/InvGauss) to ensure the most accurate and unbiased population
%  estimates.
%  ========================================================================


%% ========================================================================
%  NOTE: BROADCASTING (NETWORK)
%  ========================================================================
%  When analyzing unit-level responses (e.g., uRcv) influenced by network-
%  level properties (e.g., Dimensionality), data must be structured
%  hierarchically within the input table.
%
%  1. Implementation:
%     Each row represents a single unit. Network-level metrics are
%     duplicated for every unit belonging to that specific recording or
%     animal. The random effect term (1|Name) prevents "pseudoreplication"
%     by accounting for the non-independence of units within one culture.
%
%  2. Interaction of Scales:
%     This allows the model to test if unit-specific features (Firing Rate)
%     interact with network features (Dimensionality) to predict outcome.
%  ========================================================================

%% ========================================================================
%  NOTE: CONTINUOUS PREDICTOR SCALING
%  ========================================================================
%  Continuous predictors (e.g., fr, dimensionality) should typically be
%  standardized (Z-scored) before fitting in GLME.
%
%  1. Numerical Stability:
%     Prevents large-magnitude variables from dominating the gradient
%     descent during Maximum Likelihood Estimation.
%
%  2. Interpretability:
%     If predictors are Z-scored, the resulting coefficient (beta)
%     represents the change in log-odds per one standard deviation
%     increase in the predictor.
%
%  3. Centering:
%     Centering variables around the mean ensures that the intercept
%     represents the expected outcome at the average value of the data.
%  ========================================================================

%% ========================================================================
%  NOTE: RESPONSE VARIABLE TRANSFORMATION
%  ========================================================================
%  The decision to transform or Z-score the response variable (y) depends
%  strictly on the chosen distribution and link function.
%
%  1. Binary Response (Logistic/Binomial):
%     - DO NOT Z-score.
%     - The variable must remain binary (0 and 1).
%     - Transformation would break the Logit link function.
%
%  2. Continuous Response (Gaussian/Linear):
%     - Z-scoring is optional.
%     - Use to compare across modalities.
%     - Results in "standardized coefficients."
%
%  3. Interpretation Trade-off:
%     - Raw units: High biological meaning.
%     - Z-scored: High statistical comparability.
%
%  It is technically common to Z-score the PREDICTORS while leaving the
%  RESPONSE in its original units (or log-units).
%
%  Fixed Effects Interpretation:
%     - If x is Z-scored and y is raw:
%     - Beta = change in y-units per 1-SD change in x.
%  ========================================================================

%% ========================================================================
%  NOTE: LOG-NORMAL AND GAMMA TRANSFORMATIONS
%  ========================================================================
%  If your continuous response variable is heavily right-skewed (common in
%  firing rate distributions), a log-transform is preferred over Z-scoring.
%
%  1. Physiological Scaling:
%     - Neural firing often follows log-normal distributions.
%     - Transformation linearizes the multiplicative relationships.
%
%  2. Technical Requirement:
%     - Ensure no zero-values exist (see "Handling Zeros" note below).
%     - Log-transforming before Z-scoring is a common "best practice."
%  ========================================================================

%% ========================================================================
%  NOTE: LOG BASE SELECTION
%  ========================================================================
%  When fitting a Log-Normal model - transforming the response variable and
%  using a Gaussian distribution - the choice of logarithmic base (Base e
%  vs. Base 10) is mathematically equivalent but dictates the
%  interpretation of the fixed effects (beta).
%
%  1. Natural Log (Base e, ln): The Statistical Standard
%     In scientific research and LME modeling, the natural logarithm is the
%     nearly universal standard.
%
%     - Fractional Interpretation: For small values, the coefficient (beta)
%       represents the approximate percentage change in the response variable.
%       For example, a beta of 0.06 translates to a 6% increase in y per
%       unit increase in x.
%     - Mathematical Elegance: The derivative of ln(y) is 1/y, simplifying the
%       likelihood estimation and the recovery of geometric means via the
%       exponential function `exp()`.
%
%  2. Common Log (Base 10, log10): The Engineering Standard
%     Base 10 is preferred when data spans several orders of magnitude or when
%     the physical scale is historically logarithmic (e.g., Decibels).
%
%     - Order of Magnitude: A change of 1 on the log10 scale represents a
%       ten-fold increase in the raw units.
%     - Interpretation: Unlike natural logs, beta coefficients are not directly
%       interpretable as percentage changes; a 0.01 shift corresponds to a
%       2.3% increase, making mental approximation less intuitive.
%
%  3. Technical Selection in LME:
%     If the goal is to satisfy the assumption of normality for `fitlme`,
%     base e is the technically superior choice because it aligns with the
%     standard reporting of "percentage effects" in biological and
%     physiological literature.
%  ========================================================================





