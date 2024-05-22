
close all;
clear all;
clc;

mice_ids = [123, 126, 129, 130];
all_days = {{'221219', '221221'},...
    {'230111', '230113'},...
    {'230123', '230125'},...
    {'230322', '230324'}};

for g = 1:numel(mice_ids)

    T_prev = 0;
    mouse_id = mice_ids(g);
    mrates = [];
    days = all_days{g};
    colors = jet(numel(days));

    for ii = 1:numel(days)

        day = days{ii};
        name = sprintf('lh%d_%s',mouse_id,day);
        DATA_DIR = [sprintf('../LH_data/lh%d',mouse_id) '/' name '/'];
        load([DATA_DIR name '.spikes.cellinfo']);
        load([DATA_DIR name '.units']);

        mrate_all = [];
        slopes_all = [];
        vrate_all = [];
        dim_all = [];

        % Get firing rates
        spktms = spikes.times;

        binsize = 60; % sec
        hour = 60^2 / binsize;
        winBL = [1 Inf];
        fr = calc_fr(spktms, 'graphics', false, ...
            'binsize', binsize, 'saveVar', false, 'smet', 'MA', 'winBL', winBL);

        rates = fr.strd;
        rates = rates(units.fs',:);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        nhours = 1;
        windsz = nhours * hour; % in hours
        inc = windsz;
        T = size(rates,2);
        taxis = T_prev + [1:T]/hour;
        T_prev = taxis(end);
        figure(g); hold on;
        plot(taxis, mean(rates,1),'color',colors(ii,:)); hold on;
        xlabel('time (h)'); ylabel('MFR');
        title(mouse_id);
        set(gca,'fontsize',15);

        tottime = size(rates,2); % 24h
        nwind = floor(tottime/hour);

        Ncells = size(rates,1);
        N_all(g,ii) = Ncells;
        pvecs1 = nan(min(Ncells,inc),nwind);
        popvec = nan(Ncells,nwind);
        rdim = nan(nwind,1);
        mrate = nan(nwind,1);
        vrate = nan(nwind,1);
        dmrate = zeros(nwind,1);

        for j=1:(nwind)

            inc = windsz;
            timewindow = [(j-1) * inc + 1 : j * inc];
            Y = rates(:,timewindow);
            mrate(j) = mean(mean(Y,2));
            vrate(j) = mean(var(Y'));
            popvec(:,j) = mean(Y,2) / norm(mean(Y,2));
        end

        PVcorr = corr(popvec);
        nvecs = size(pvecs1,2);
        mv = nan(nvecs-1,1);

        figure(1010*g + 100*ii); hold on;
        for kk = 1:nvecs-1

            v = diag(PVcorr,kk);
            mv(kk) = mean(v);
            plot(kk*nhours,mv(kk),'.k','MarkerSize',15); hold on;
            plot(kk*nhours,v,'.','color',[.5 .5 .5]);

        end
        set(gca,'FontSize',15); axis square;
        xlabel('time elapsed [h]'); ylabel('PV correlation');

        tmv = mv';
        mv_all{ii} = tmv;

        xaxis = [1:nvecs-1]*nhours;
        p = polyfit(xaxis,mv,1);
        plot(xaxis,p(2) + p(1)*xaxis,'b');
        title(sprintf('slope = %.3f',p(1)));
        drift_rates(g,ii) = p(1);
        plot_nwind = length(mv);

        figure(1000*g); title(mouse_id);
        plot(xaxis(1:plot_nwind-1),mv(1:plot_nwind-1),'LineWidth',2,'color',colors(ii,:)); hold on;
        set(gca,'FontSize',15); axis square; box off;
        xlabel('time elapsed [h]'); ylabel('PV correlation');
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

figure;
for g=1:numel(days)

    plot(g,abs(drift_rates(:,g)),'.','markersize',15,'color',colors(g,:)); hold on; box off;
    aa = abs(drift_rates(:,g));
    m = mean(aa);
    plot(g,m,'.k','MarkerSize',35);

end

xlim([0 numel(days)+1]);
set(gca,'FontSize',15);
ylabel('drift rate [1/h]');
xticks([1:numel(days)]);
xticklabels({'Baseline', 'MK-801 day 2'});
