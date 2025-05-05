function [r_mean, r_sem, fh] = mcu_FRvBL(lme_tbl)
% analyzes correlation between bout length and firing rate per neuron
%
% INPUT
%   lme_tbl     table with fields FR, BoutDur, Group, State, UnitID
%
% OUTPUT
%   r_mean      mean correlation coefficient [state x group]
%   r_sem       standard error of correlation [state x group]
%   fh          figure handle showing distribution of correlations
%
% 08 Jan 24

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get unique states, groups, and units
states = unique(lme_tbl.State);
groups = double(unique(lme_tbl.Group));
units = unique(lme_tbl.UnitID);
nstates = length(states);
ngroups = length(groups);
nunits = length(units);

% initialize unit correlations
r_units = nan(nunits, nstates);   % correlation coef per unit
grp_units = nan(nunits, 1);       % group assignment per unit

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate correlations per unit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iunit = 1 : nunits
    % get all data from current unit
    unit_tbl = lme_tbl(lme_tbl.UnitID == units(iunit), :);
    
    % store group
    grp_units(iunit) = double(unit_tbl.Group(1));
    
    % calculate correlation per state
    for istate = 1 : nstates
        % get state data
        idx = unit_tbl.State == states(istate);
        if sum(idx) > 5     % minimum 5 observations
            x = unit_tbl.BoutDur(idx);
            y = log(unit_tbl.FR(idx));
            
            % remove inf/nan
            validIdx = isfinite(x) & isfinite(y);
            if sum(validIdx) > 5
                [r, ~] = corr(x(validIdx), y(validIdx));
                r_units(iunit, istate) = r;
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate stats
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize output
r_mean = nan(nstates, ngroups);
r_sem = nan(nstates, ngroups);

% stats per state and group
for istate = 1 : nstates
    for igroup = 1 : ngroups
        idx = (grp_units) == groups(igroup);
        r_mean(istate, igroup) = nanmean(r_units(idx, istate));
        r_sem(istate, igroup) = nanstd(r_units(idx, istate)) / ...
            sqrt(sum(~isnan(r_units(idx, istate))));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% visualization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fh = figure;
colors = {'b', 'r'};
for istate = 1 : nstates
    subplot(1, nstates, istate)
    hold on
    
    % violin plot per group
    for igroup = 1 : ngroups
        idx = grp_units == groups(igroup);
        r_grp = r_units(idx, istate);
        
        % remove nan for plotting
        r_grp = r_grp(~isnan(r_grp));
        
        % plot distribution
        violin(r_grp, igroup, 'FaceColor', colors{igroup}, 'EdgeColor', colors{igroup});
        
        % add mean and sem
        errorbar(igroup, r_mean(istate, igroup), r_sem(istate, igroup),...
            'k', 'LineWidth', 2);
    end
    
    % aesthetics
    xlabel('Group')
    ylabel('Correlation (r)')
    title(sprintf('State %d', istate))
    box off
    ylim([-1 1])
end


end

% local function for violin plot if not already available
function violin(data, pos, varargin)
    % very basic violin plot
    [f, xi] = ksdensity(data);
    f = f / max(f) * 0.3;    % scale width
    
    % plot violin
    fill([f; -f] + pos, [xi; flipud(xi)], 'b', varargin{:});
    
    % add median line
    line([-0.3 0.3] + pos, [median(data) median(data)], 'Color', 'k', 'LineWidth', 2);
end


