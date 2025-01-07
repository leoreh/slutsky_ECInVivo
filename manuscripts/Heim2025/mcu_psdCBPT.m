function [stats, cfg_cbp] = mcu_psdCBPT(psd_data, psd_cfg)

% wrapper for fieldtrip's ft_freqstatistics for analyzing psd data from
% mcu_psdOrg
%
% INPUT
%   psd_data      cell of psd mats organized by mcu_psdOrg
%   psd_cfg       config from mcu_psdOrg
%
% OUTPUT
%   stats     output from ft_freqstatistics
%   cfg_out   config used for ft_freqstatistics
%
% CALLS
%   ft_freqstatistics
%
% 05 Jan 25 LH


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preperations 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% exclude pli frequencies
flg_pli = true;
faxis = psd_cfg.faxis;
if flg_pli
    freq_mask = ~(faxis >= 48 & faxis <= 52);  % Logical mask for frequencies to keep
else
    freq_mask = true(length(faxis), 1);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% organize data for fieldtrip based on experimental design
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% experiment specific parameters
switch psd_cfg.expDsgn
    case 1                              % '{state}[mouse x freq x time]'
        dimord = 'subj_chan_freq_time';
        
        % dependent samples (same mouse different states)
        cfg_cbp.statistic = 'ft_statfun_depsamplesT';
        design = zeros(2, 2*size(psd_data{1}, 1));
        design(1, :) = [1 : size(psd_data{1}, 1) 1 : size(psd_data{1}, 1)];
        design(2, :) = [ones(1, size(psd_data{1}, 1)) 2*ones(1, size(psd_data{1}, 1))];
        cfg_cbp.uvar = 1;       % unit variable (subject number)
        cfg_cbp.ivar = 2;       % independent variable (state)

    case 2                              % '{gen}[freq x mouse]'        
        dimord = 'subj_chan_freq';

        % independent samples (different mice)
        cfg_cbp.statistic = 'ft_statfun_indepsamplesT';
        design = [ones(1, size(psd_data{1}, 2)),...
            2 * ones(1, size(psd_data{2}, 2))];
        cfg_cbp.ivar = 1;

    otherwise
        error('Unsupported data format')
end

% data structure (per group)
for igrp = 1 : 2
    
    switch psd_cfg.expDsgn
        case 1                      % [mouse chan day freq]
            psdMat = psd_data{igrp}(:, freq_mask, :);
            powspctrm = permute(psdMat, [1 4 2 3]);     
            freq_data(igrp).time = 1 : size(psd_data{1}, 3);

        case 2                      % [subj chan freq]
            psdMat = psd_data{igrp}(freq_mask, :);
            powspctrm = permute(psdMat', [1 3 2]);      
    end
    freq_data(igrp).powspctrm = powspctrm;
    freq_data(igrp).freq = faxis(freq_mask);
    freq_data(igrp).dimord = dimord;
    freq_data(igrp).label = {'dummy'};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% configure statistics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cfg_cbp.method = 'montecarlo';
cfg_cbp.correctm = 'tfce';
cfg_cbp.clusteralpha = 0.05;
cfg_cbp.clusterstatistic = 'maxsum';
cfg_cbp.minnbchan = 0;
cfg_cbp.tail = 0;
cfg_cbp.alpha = 0.05;
cfg_cbp.numrandomization = 1000;
cfg_cbp.design = design;
cfg_cbp.neighbours = []; % no spatial clustering since we only have one channel

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run statistics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stats = ft_freqstatistics(cfg_cbp, freq_data(1), freq_data(2));

end

% EOF

% performs a cluster-based permutation test to compare PSDs between two
% groups using FieldTrip. Cluster-based permutation testing is a
% non-parametric approach designed to address the multiple comparisons
% problem that arises in electrophysiological data, such as EEG, MEG, and
% LFP signals. When comparing PSDs across many frequency bins, performing
% independent t-tests at each bin can lead to inflated false positives.
% Traditional correction methods, like Bonferroni, are often overly
% conservative, reducing sensitivity to real effects. Cluster-based
% permutation circumvents this by leveraging the spatial or spectral
% structure of the data, operating under the assumption that meaningful
% effects are likely to occur across contiguous bins rather than isolated
% points.

% The process begins by performing pointwise t-tests (or other statistical
% tests) at each frequency bin, comparing PSDs between the two groups.
% Frequency bins that exceed a predefined threshold (e.g., p < 0.05) are
% marked as significant. Adjacent bins showing significance are grouped
% into clusters, and each cluster is assigned a summary statistic,
% typically the sum of the t-values across all bins in the cluster. This
% summary value represents the strength and extent of the difference
% between groups across contiguous frequencies. 

% To assess the statistical significance of the observed clusters, the data
% labels (i.e., group memberships) are randomly permuted many times (e.g.,
% 1000 permutations) to generate a null distribution of cluster statistics.
% In each permutation, the same t-tests and clustering procedures are
% applied, yielding a distribution of the largest cluster statistics
% expected by chance. The observed clusters from the original (unpermuted)
% data are then compared to this null distribution. A cluster is deemed
% statistically significant if its size or strength exceeds 95% of the
% clusters from the permuted data (corresponding to p < 0.05).

% A notable advantage of this approach is its ability to preserve
% sensitivity to true effects while maintaining strict control over false
% positives. By testing clusters rather than individual bins, cluster-based
% permutation reduces the number of comparisons without discarding
% meaningful information. The method is particularly effective for
% detecting differences in oscillatory power, evoked potentials, or
% time-frequency responses, making it a standard tool in neurophysiology
% and cognitive neuroscience. The FieldTrip implementation of this method,
% as used in the provided function, adheres to the principles laid out by
% Maris and Oostenveld (2007), a seminal paper that established the
% framework for cluster-based permutation testing in EEG and MEG data.

% Another option for when running this test is to set cfg.correctm =
% 'tfce'. Threshold-Free Cluster Enhancement (TFCE) is a statistical method
% developed to address limitations in traditional cluster-based statistics.
% In conventional cluster-based approaches, you need to set an arbitrary
% threshold to form clusters - any data points above this threshold are
% grouped into clusters, while those below are discarded. This binary
% decision can be problematic because it makes the results heavily
% dependent on your chosen threshold. A slightly different threshold might
% give very different results, and potentially meaningful sub-threshold
% effects might be completely missed.

% TFCE solves this by eliminating the
% need for this arbitrary threshold. Instead, it considers every possible
% threshold and combines information about both the height (strength) of
% the statistical values and their spatial extent (how many neighboring
% points show similar effects). For each data point, TFCE calculates a
% score by integrating across all possible thresholds below its height,
% weighting both how far above the threshold the data point is (height) and
% how large the cluster would be at that threshold (extent). This
% integration means that a data point gets support from both strong
% neighboring effects and weaker but more extensive effects.

% The method has two main parameters: H (height) and E (extent) exponents
% that control how these factors are weighted. The default values (H=2,
% E=0.5) have been found to work well across many applications, but you can
% adjust them if you have specific reasons to weight either height or
% extent differently. Higher H values give more weight to intense effects,
% while higher E values emphasize spatial extent. 
% 
% The trade-off is increased computational time, as TFCE needs to calculate
% scores across many thresholds, but this is usually worth it for the
% increased sensitivity and robustness of the results.
%
% references
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% https://www.fieldtriptoolbox.org/tutorial/cluster_permutation_freq/#introduction
% https://www.fieldtriptoolbox.org/faq/clusterstats_interpretation/
% Maris E., Oostenveld R. Nonparametric statistical testing of EEG- and MEG-data. J Neurosci Methods. 2007 Apr 10;