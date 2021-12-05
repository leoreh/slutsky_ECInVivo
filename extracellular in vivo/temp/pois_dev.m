

st_metrics = ["doublets"; "royer"; "lidor"; "mizuseki"; "cv"; "cv2"; "lv"; "lvr"];
nmetrics = length(st_metrics);
nsub = numSubplots(nmetrics);
frBoundries = [0 Inf; 0 Inf];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load and arrange data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% relavent sessions
basepaths{1} = 'D:\Data\lh93\lh93_210811_102035';
basepaths{2} = 'K:\Data\lh95\lh95_210826_080100';
basepaths{3} = 'G:\lh87\lh87_210523_100607';
basepaths{4} = 'G:\lh81\lh81_210206_044300';

% load
varArray = getSessionVars('dirnames', basepaths, 'mousepath', [],...
    'sortDir', false);

% initialize
for imetric = 1 : nmetrics
    st_cat.(st_metrics{imetric}) = [];
end
pyrUnits = logical([]); intUnits = logical([]);
fr_wake = []; fr_nrem = []; stimes = []; mfr = []; 
acg_wide = []; acg_narrow = [];

% cat
for isession = 1 : length(basepaths)
    assignVars(varArray, isession)
    for imetric = 1 : nmetrics
        st_cat.(st_metrics{imetric}) = [st_cat.(st_metrics{imetric});...
            st.(st_metrics{imetric})];
    end
    acg_wide = [acg_wide, st.acg_wide(:, :, 1)];
    acg_narrow = [acg_narrow, st.acg_narrow(:, :, 1)];
    pyrUnits = [pyrUnits; selectUnits(spikes, cm, fr, 1, [1 : 4], frBoundries, 'pyr')];
    intUnits = [intUnits; selectUnits(spikes, cm, fr, 1, [1 : 4], frBoundries, 'int')];
    fr_wake = [fr_wake; mean(fr.states.fr{1}, 2)];
    fr_nrem = [fr_nrem; mean(fr.states.fr{4}, 2)];
    mfr = [mfr; fr.mfr];
    stimes = [stimes, spikes.times];
end
stimes = stimes(:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TO DO
% movsum instead of histcounts on spktimes 
% diff and params of obsv_dist

nunits = length(stimes);
sunits = 1 : nunits;
sunits = 10 : 15;
graphics = true;

for iunit = sunits
    
    % calc
    expc_dist = makedist('Poisson', 'lambda', mfr(iunit));
%     spktimes_edges = [0 : mfr(iunit) : max_len];
%     spkcounts_mfr = histcounts(stimes{iunit}, spktimes_edges);
    spkcounts_mfr = histcounts(stimes{iunit}, [0 : max(stimes{iunit})]);
    [obsv_dist, pois_edges] = histcounts(spkcounts_mfr, 49, 'Normalization', 'probability');
    x_val = [0 : 1 : pois_edges(end)];
    obsv_dist = interp1(pois_edges(1 : end - 1), obsv_dist, x_val, 'spline', 'extrap');
    obsv_dist = obsv_dist / sum(obsv_dist);
    y_val = pdf(expc_dist, x_val);
    y_val = y_val / sum(y_val);
    
    pois_dev(iunit) = sum((obsv_dist - y_val).^2);
    
    % Kullback–Leibler divergence
    kl(iunit) = abs(kldiv(x_val, obsv_dist + eps, y_val + eps));
    
    % fano factor
    ff(iunit) = var(spkcounts_mfr) / mean(spkcounts_mfr);
    
    % fig per unit --------------------------------------------------------
    if graphics
        fh = figure('Visible', 'on');
        subplot(2, 2, 1)
        plotCCG('ccg', acg_narrow(:, iunit),...
            't', st.info.acg_narrow_tstamps);
        ylabel('Rate [Hz]')
        xlabel('Time [ms]')
        subtitle(sprintf('fr: %.2f', mfr(iunit)))
        
        subplot(2, 2, 2)
        plotCCG('ccg', acg_wide(:, iunit),...
            't', st.info.acg_wide_tstamps);
        ylabel('Rate [Hz]')
        xlabel('Time [ms]')
        subtitle(sprintf('burstiness: %.2f', st_cat.lidor(iunit, 1)))
        
        subplot(2, 2, 3)
        x_max = 0.1;
        isi = diff(stimes{iunit});
        histogram(isi(isi < x_max), 100, 'Normalization', 'probability')
        xlim([0 x_max])
        xlabel('Time [ms]')
        subtitle(sprintf('lvr: %.2f', st_cat.lvr(iunit, 1)))
        
        subplot(2, 2, 4)
        plot(x_val, y_val)
        hold on
        plot(x_val, obsv_dist)
        legend({'Expc', 'Obsv'})
        xlabel('Spks per bin')
        subtitle(sprintf('dev=%.2f; kl=%.2f; ff=%.2f',...
            pois_dev(iunit), kl(iunit), ff(iunit)))
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% corr w/ fr

fh = figure;

sb1 = subplot(2, 2, 1);
x1 = fr_wake(pyrUnits);
x2 = fr_wake(intUnits);
y1 = pois_dev(pyrUnits)';
y2 = pois_dev(intUnits)';
corrFR(x1, y1, x2, y2, sb1)
title('Diff^2')

sb2 = subplot(2, 2, 2);
x1 = fr_wake(pyrUnits);
x2 = fr_wake(intUnits);
y1 = kl(pyrUnits)';
y2 = kl(intUnits)';
corrFR(x1, y1, x2, y2, sb2)
title('Kullback–Leibler')

sb3 = subplot(2, 2, 3);
x1 = fr_wake(pyrUnits);
x2 = fr_wake(intUnits);
y1 = ff(pyrUnits)';
y2 = ff(intUnits)';
corrFR(x1, y1, x2, y2, sb3)
title('Fano Factor')

% corr with other metrics -------------------------------------------------

fh = figure;
sz = log10(fr_wake) * 30;
sz(sz < 8) = 8;

subplot(2, 2, 1)
pyr_x = st_cat.lvr(pyrUnits, 1); % wake
pyr_y = kl(pyrUnits); % nrem
int_x = st_cat.lvr(intUnits, 1);
int_y = kl(intUnits);
scatter(pyr_x, pyr_y, sz(pyrUnits), 'b', 'filled')
hold on
scatter(int_x, int_y, sz(intUnits), 'r', 'filled')
xlabel('LvR')
ylabel('KL')

subplot(2, 2, 2)
pyr_x = st_cat.lidor(pyrUnits, 1); % wake
pyr_y = kl(pyrUnits); % nrem
int_x = st_cat.lidor(intUnits, 1);
int_y = kl(intUnits);
scatter(pyr_x, pyr_y, sz(pyrUnits), 'b', 'filled')
hold on
scatter(int_x, int_y, sz(intUnits), 'r', 'filled')
xlabel('Burstiness')
ylabel('KL')

subplot(2, 2, 3)
pyr_x = pois_dev(pyrUnits); % wake
pyr_y = kl(pyrUnits); % nrem
int_x = pois_dev(intUnits);
int_y = kl(intUnits);
scatter(pyr_x, pyr_y, sz(pyrUnits), 'b', 'filled')
hold on
scatter(int_x, int_y, sz(intUnits), 'r', 'filled')
xlabel('Pois Dev')
ylabel('KL')

subplot(2, 2, 4)
pyr_x = st_cat.lidor(pyrUnits, 1); % wake
pyr_y = st_cat.lvr(pyrUnits); % nrem
int_x = st_cat.lidor(intUnits, 1);
int_y = st_cat.lvr(intUnits);
scatter(pyr_x, pyr_y, sz(pyrUnits), 'b', 'filled')
hold on
scatter(int_x, int_y, sz(intUnits), 'r', 'filled')
xlabel('Burstiness')
ylabel('LVR')