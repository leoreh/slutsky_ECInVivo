
% analyze anesthesia states. Is essentialy a wrapper for analysis at the
% mouse (aneStates_m) and group (aneStates_g) level.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
basepath{1} = 'E:\Data\Others\DZ\Field\Data\APPPS1_short';
basepath{2} = 'E:\Data\Others\DZ\Field\Data\APPPS1_long';
basepath{3} = 'E:\Data\Others\DZ\Field\Data\WT_long';
basepath{4} = 'E:\Data\Others\DZ\Field\Data\WT_short';
rm = cell(4, 1);

forceA = true;
forceL = false;
saveFig = true;
graphics = true;
saveVar = true;
ch = 1;
smf = 6;
fs = 1250;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% group data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1 : length(basepath)
    as{i} = aneStates_g('basepath', basepath{i}, 'rm', rm{i},...
        'graphics', graphics, 'saveVar', saveVar,...
        'saveFig', saveFig, 'forceA', forceA);
    
    nspks{i} = as{i}.nspks;
    nspksEp{i} = as{i}.nspksEp;
    recDur{i} = as{i}.recDur;
    epDur{i} = as{i}.epDur;
    thr{i} = as{i}.thr;
    grp{i} = as{i}.grp;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c = [1 0 0; 1 0 1; 0 0 1; 0 1 1];
c2 = 'rmbc';
if graphics
    
    figure
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    
    % duration bar
    subplot(3, 3, 1 : 3)
    hold on
    for i = [1, 4]
        stdshade(as{i}.iis, 0.5, c2(i), as{i}.t / fs / 60)
    end
    axis tight
    ylabel('Rate [IIS / min]')
    xlim([1 50])
    box off
    set(gca, 'TickLength', [0 0])
    title('Short Recordings')
    
    % duration bar
    subplot(3, 3, 4 : 6)
    hold on
    for i = [2, 3]
        stdshade(as{i}.iis, 0.5, c2(i), as{i}.t / fs / 60)
    end
    axis tight
    ylabel('Rate [IIS / min]')
    xlabel('Time [m]')
    xlim([1 100])
    box off
    set(gca, 'TickLength', [0 0])
    title('Long Recordings')
    
    % duration bar
    subplot(3, 3, 7)
    b = bar(([cellfun(@mean, epDur)' (cellfun(@mean, recDur) - cellfun(@mean, epDur))']), 'stacked',...
        'FaceColor', 'flat');
    for i = 1 : 4
        b(1).CData(i, :) = c(i, :) + 0.6;
        b(2).CData(i, :) = c(i, :);
    end
    axis tight
    xticklabels(grp)
    legend({'BS', 'Not BS'})
    ylabel('Duration [m]')
    title('Duration')
    box off
    set(gca, 'TickLength', [0 0])
    
    % IIS bar
    subplot(3, 3, 8)
    b = bar(([cellfun(@mean, nspksEp)' (cellfun(@mean, nspks) - cellfun(@mean, nspksEp))']), 'stacked',...
        'FaceColor', 'flat');
    for i = 1 : 4
        b(1).CData(i, :) = c(i, :) + 0.6;
        b(2).CData(i, :) = c(i, :);
    end
    axis tight
    xticklabels(grp)
    legend({'IIS in BS', 'IIS out BS'})
    ylabel('IIS [#]')
    title('IIS')
    box off
    set(gca, 'TickLength', [0 0])
    
    % threshold bar
    subplot(3, 3, 9)
    b = bar(cellfun(@mean, thr)', 'FaceColor', 'flat');
    for i = 1 : 4
        b(1).CData(i, :) = c(i, :);
    end
    axis tight
    xticklabels(grp)
    ylabel('Threshold [mV]')
    title('Threshold')
    box off
    set(gca, 'TickLength', [0 0])
    
    if saveFig
        figname = ['summary'];
        export_fig(figname, '-tif', '-transparent')
        % savePdf(figname, basepath, ff)
    end
end
