function stat = stat_cluPermutation(dataCell, faxis)

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

% Another option for when running this test is to set cfg.correctm = 'tfce';
% Threshold-Free Cluster Enhancement (TFCE) is a statistical method that
% improves sensitivity in cluster-based analyses by eliminating the need to
% predefine a cluster-forming threshold. Traditional cluster-based
% permutation testing relies on setting a significance threshold (e.g., p <
% 0.05) to identify contiguous bins for clustering. However, this approach
% can miss smaller but meaningful effects if they do not exceed the
% arbitrary threshold. TFCE addresses this by enhancing signal intensity
% across all bins based on both their magnitude (t-value) and
% spatial/spectral extent (cluster size), without requiring a hard cutoff.
% Bins with stronger effects and larger cluster sizes receive greater
% enhancement, while isolated bins with weak effects receive minimal
% enhancement. This results in a continuous, weighted representation of
% cluster-like structures, increasing the likelihood of detecting true
% effects without inflating false positives. TFCE has become widely used in
% neuroimaging and electrophysiology for analyzing time-frequency data,
% power spectra, and spatial EEG/MEG signals, providing a more robust
% alternative to standard cluster-based correction methods.
% 
% The input dataCell is a cell array containing PSD data matrices for each
% group. Each matrix is organized as [frequency x subjects], where the rows
% represent frequency bins and the columns correspond to individual
% subjects or trials. The second input, faxis, is a vector specifying the
% frequency axis (e.g., 1 to 652 Hz). This vector aligns with the frequency
% dimension of the PSD data, ensuring correct mapping between statistical
% results and frequency bins.
% 
% The function first sets up the FieldTrip configuration (cfg) for a Monte
% Carlo permutation test. The method is specified as 'montecarlo' and the
% test statistic as 'indepsamplesT', indicating the use of an independent
% samples t-test, suitable for comparing two groups of unpaired data. The
% cluster-based correction method ('cluster') is applied to control for
% multiple comparisons, and a threshold of 0.05 (cfg.clusteralpha) is used
% to determine initial cluster formation. The function is two-tailed
% (cfg.tail = 0), allowing it to detect both positive and negative
% differences between groups. The number of random permutations is set to
% 1000, providing sufficient iterations to build a robust null
% distribution.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% references 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% https://www.fieldtriptoolbox.org/tutorial/cluster_permutation_freq/#introduction
% https://www.fieldtriptoolbox.org/faq/clusterstats_interpretation/
% Maris E., Oostenveld R. Nonparametric statistical testing of EEG- and MEG-data. J Neurosci Methods. 2007 Apr 10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preperations 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ngrps = length(dataCell);
ndims_cell = unique(cellfun(@(x) ndims(x), dataCell, 'uni', true));

% exclude pli frequencies
flg_pli = true;
if flg_pli
    freq_mask = ~(faxis >= 48 & faxis <= 52);  % Logical mask for frequencies to keep
else
    freq_mask = true(length(faxis), 1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% study configuration 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cfg = [];
cfg.method = 'montecarlo';
cfg.statistic = 'indepsamplesT';
cfg.correctm = 'cluster';
cfg.clusteralpha = 0.05;
cfg.tail = 0;
cfg.clustertail = 0;
cfg.numrandomization = 1000;
cfg.neighbours = [];  % No spatial neighbors

% study design
if ndims_cell == 2
    nmice = cellfun(@(x) size(x, 2), dataCell, 'uni', true);
    cfg.design = [ones(1, nmice(1)), 2 * ones(1, nmice(2))];  % group labels
    cfg.ivar = 1;  % Independent variable (state)

elseif ndims_cell == 3

    cfg.statistic = 'depsamplesT';
    nmice = unique(cellfun(@(x) size(x, 3), dataCell, 'uni', true));
    ndays = unique(cellfun(@(x) size(x, 2), dataCell, 'uni', true));
    cfg.design = [1:nmice, 1:nmice;             % Mouse ID (uvar - unit of observation)
        ones(1, nmice), 2 * ones(1, nmice)];    % State (ivar - condition)

    cfg.uvar = 1;  % Unit of observation (mouse)
    cfg.ivar = 2;  % Independent variable (state)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% organize data structure 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create ft structure
clear freq_strct
for igrp = 1 : ngrps
    

    if ndims_cell == 2
        dataMat = dataCell{igrp}(freq_mask, :);
        freq_strct(igrp).powspctrm(:, 1, :) = dataMat';
        freq_strct(igrp).dimord = 'subj_freq';

    elseif ndims_cell == 3
        dataMat = dataCell{igrp}(freq_mask, :, :);
        tempData = permute(dataMat, [3, 1, 2]);
        freq_strct(igrp).powspctrm = reshape(tempData, size(tempData, 1), 1, size(tempData, 2), size(tempData, 3));

        freq_strct(igrp).dimord = 'subj_chan_freq_time';
        freq_strct(igrp).time = 1 : ndays;
    end

    freq_strct(igrp).freq = faxis(freq_mask);
    freq_strct(igrp).label = {'dummy'};
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run statistical comparison
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stat = ft_freqstatistics(cfg, freq_strct(1), freq_strct(2));

end
