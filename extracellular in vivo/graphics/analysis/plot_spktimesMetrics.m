
% the following assumes st was computed for awake and sleep

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
fr_wake = []; fr_nrem = []; stimes = [];

% cat
for isession = 1 : length(basepaths)
    assignVars(varArray, isession)
    for imetric = 1 : nmetrics
        st_cat.(st_metrics{imetric}) = [st_cat.(st_metrics{imetric});...
            st.(st_metrics{imetric})];
    end
    pyrUnits = [pyrUnits; selectUnits(spikes, cm, fr, 1, [1 : 4], frBoundries, 'pyr')];
    intUnits = [intUnits; selectUnits(spikes, cm, fr, 1, [1 : 4], frBoundries, 'int')];
    fr_wake = [fr_wake; mean(fr.states.fr{1}, 2)];
    fr_nrem = [fr_nrem; mean(fr.states.fr{4}, 2)];
    stimes = [stimes, spikes.times];
end
stimes = stimes(:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% corr w/ fr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fh = figure;
xLimit = [5 * 1e-2 5 * 1e1];
for imetric = 1 : nmetrics
    subplot(nsub(1), nsub(2), imetric)
    
    % data
    pyr_x = fr_wake(pyrUnits);
    pyr_y = st_cat.(st_metrics{imetric})(pyrUnits, 1);
    int_x = fr_wake(intUnits);
    int_y = st_cat.(st_metrics{imetric})(intUnits, 1);    
    if any(strcmp(st_metrics{imetric}, {'royer', 'mizuseki', 'cv'}))
        pyr_y = log10(pyr_y);
        int_y = log10(int_y);
    end
    
    % corr and linfit
    [pyr_r, pyr_p] = corr(log10(pyr_x), pyr_y);
    [int_r, int_p] = corr(log10(int_x), int_y);
    pyr_lf = polyfit(log10(pyr_x), pyr_y, 1);
    pyr_y2 = polyval(pyr_lf, log10(pyr_x));
    int_lf = polyfit(log10(int_x), int_y, 1);
    int_y2 = polyval(int_lf, log10(int_x));
    
    % plot
    scatter(pyr_x, pyr_y, 20, 'b', 'filled')
    hold on
    scatter(int_x, int_y, 20, 'r', 'filled')
    plot(pyr_x, pyr_y2, '--b')
    plot(int_x, int_y2, '--r')
    set(gca, 'XScale', 'log')
    hold on
    xlim(xLimit)
    xticks([1e-2, 1e-1, 1e1])
    xlabel('Mean Firing Rate [Hz]')
    title(st_metrics{imetric})
    subtitle(sprintf('pyr: r=%.2f, p=%.2f \nint: r=%.2f, p=%.2f',...
        pyr_r, pyr_p, int_r, int_p))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WAKE vs NREM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fh = figure;
for imetric = 1 : nmetrics
    subplot(nsub(1), nsub(2), imetric)
    pyr_x = st_cat.(st_metrics{imetric})(pyrUnits, 1); % wake
    pyr_y = st_cat.(st_metrics{imetric})(pyrUnits, 2); % nrem
    int_x = st_cat.(st_metrics{imetric})(intUnits, 1); 
    int_y = st_cat.(st_metrics{imetric})(intUnits, 2); 
    sz = log10(fr_wake) * 30;
    sz(sz < 8) = 8;
    scatter(pyr_x, pyr_y, sz(pyrUnits), 'b', 'filled')
    hold on
    scatter(int_x, int_y, sz(intUnits), 'r', 'filled')  
    if any(strcmp(st_metrics{imetric}, {'royer', 'mizuseki', 'doublets'}))
        set(gca, 'YScale', 'log', 'XScale', 'log')
    end
    hold on
    axLimit = [min([xlim, ylim]) max([xlim, ylim])];
    xlim(axLimit)
    ylim(axLimit)
    plot([axLimit(1), axLimit(2)], [axLimit(1), axLimit(2)], '--k')
    xlabel('WAKE')
    ylabel('NREM')
    title(st_metrics{imetric})
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% acg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fh = figure;
iunit = 3;
nwin = size(st.acg_narrow, 3);

subplot(1, nwin, 1)
plotCCG('ccg', squeeze(st.acg_narrow(:, iunit, 1)),...
    't', st.info.acg_narrow_tstamps)
ylabel('Rate [Hz]')
xlabel('Time [ms]')
title('WAKE')

subplot(1, nwin, 2)
plotCCG('ccg', squeeze(st.acg_narrow(:, iunit, 2)),...
    't', st.info.acg_narrow_tstamps)
ylabel('Rate [Hz]')
xlabel('Time [ms]')
title('NREM')
