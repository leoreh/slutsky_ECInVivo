function [idx_eq, fh] = mcu_eqBLen(lme_tbl)
% matches Bout Dur distributions between groups by selecting equal number of
% bouts per bin from each group. Uses log-spaced bins to account for the
% skewed distribution of bout lengths.
%
% Strategy:
% 1. For each state:
%    - Create log-spaced bins covering the range of bout lengths
%    - Count bouts per bin for each group
%    - For each bin, find minimum count across groups
%    - Randomly select that many bouts from each group
% 2. This ensures:
%    - Equal representation across the range despite skewed distribution
%    - Matching Bout Dur distributions between groups
%    - Random selection within each bin
%
% INPUT
%   lme_tbl     table with fields BoutDur, Group, State from mcu_lmeOrg
%
% OUTPUT
%   idx_eq      logical vector size of lme_tbl for equal bout lengths  
%   fh          figure handle comparing original and matched distributions
%
% 08 Jan 24

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get unique states and groups
states = unique(lme_tbl.State);
groups = unique(lme_tbl.Group);
nstates = length(states);
ngroups = length(groups);

% initialize output
idx_eq = false(height(lme_tbl), 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% match distributions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% settings
nbins = 50;

% process each state
for istate = 1 : nstates
   % get data for current state
   state_tbl = lme_tbl(lme_tbl.State == states(istate), :);
   
   % get bout lengths per group
   bDur = cell(ngroups, 1);
   for igroup = 1 : ngroups
       bDur{igroup} = state_tbl.BoutDur(state_tbl.Group == groups(igroup));
   end
   
   % create log-spaced edges
   minval = min(vertcat(bDur{:}));
   maxval = max(vertcat(bDur{:}));
   edges = exp(linspace(log(minval), log(maxval), nbins + 1));
   
   % count bouts per bin per group
   counts = zeros(ngroups, nbins);
   for igroup = 1 : ngroups
       counts(igroup, :) = histcounts(bDur{igroup}, edges, 'Normalization', 'count');
   end
   
   % find minimum count per bin
   min_count = min(counts, [], 1);
   
   % select bouts
   for igroup = 1 : ngroups
       % get group data
       grp_idx = state_tbl.Group == groups(igroup);
       curr_blen = state_tbl.BoutDur(grp_idx);
       
       % initialize selection
       sel_idx = false(size(curr_blen));
       
       % select per bin
       for ibin = 1 : nbins
           bin_idx = curr_blen >= edges(ibin) & curr_blen < edges(ibin + 1);
           if sum(bin_idx) > 0
               % randomly select required number of bouts
               curr_idx = find(bin_idx);
               rand_idx = randperm(length(curr_idx), min(min_count(ibin), length(curr_idx)));
               sel_idx(curr_idx(rand_idx)) = true;
           end
       end
       
       % update output
       idx_eq(lme_tbl.State == states(istate) & lme_tbl.Group == groups(igroup)) = sel_idx;
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% visualization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create figure
fh = figure;
colors = {'b', 'r'};

% process each state
for istate = 1 : nstates
   % get original distributions
   state_tbl = lme_tbl(lme_tbl.State == states(istate), :);
   bDur = cell(ngroups, 1);
   for igroup = 1 : ngroups
       bDur{igroup} = state_tbl.BoutDur(state_tbl.Group == groups(igroup));
   end
   
   % plot original distributions
   subplot(2, nstates, istate)
   hold on
   for igroup = 1 : ngroups
       histogram(bDur{igroup}, edges, 'Normalization', 'count',...
           'FaceColor', colors{igroup}, 'FaceAlpha', 0.3)
   end
   set(gca, 'XScale', 'log')
   title(sprintf('State %d - Original', istate))
   xlabel('Bout Dur (log scale)')
   ylabel('Probability')
   legend(cellstr(groups))
   
   % plot matched distributions
   subplot(2, nstates, istate + nstates)
   hold on
   for igroup = 1 : ngroups
       curr_idx = lme_tbl.State == states(istate) & ...
           lme_tbl.Group == groups(igroup) & idx_eq;
       histogram(lme_tbl.BoutDur(curr_idx), edges, 'Normalization', 'count',...
           'FaceColor', colors{igroup}, 'FaceAlpha', 0.3)
   end
   set(gca, 'XScale', 'log')
   title(sprintf('State %d - Matched', istate))
   xlabel('Bout Dur (log scale)')
   ylabel('Probability')
   legend(cellstr(groups))
end

end