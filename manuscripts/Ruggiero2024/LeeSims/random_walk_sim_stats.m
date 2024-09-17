
%  This script simulates two different systems of random processes:
%  1) Processes are independent
%  2) Processes are strongly correlated
%  
%  and computes the dimensionality of the dynamics, and the stability of
%  the population-averaged activity (Suplementary Figure 1E-H).


close all; clear all; clc;

tic

% network and simulation parameters
N = 150;
T = 200;
dt = .1;
S = 2;
Nr = 500;

dim_all = nan(S,Nr);
stability_all = nan(S,Nr);

for r = 1:Nr
    dim_r = nan(S,1);
    stability_r = nan(S,1);

    for ss = 1:S

        Tt = T/dt;
        transient = Tt/4;
        X = nan(N,Tt);
        X(:,1) = rand(N,1);
        v = ones(N,1) + randn(N,1)/sqrt(N);

        % simulate via Euler forward integration
        for i=1:Tt-1

            x = X(:,i);
            if ss == 2  % strongly correlated units
                dx = v*randn(1,1);
            else        % independent units
                dx = randn(N,1);
            end
            X(:,i+1) = x + dx*dt;

        end

        % scale and add white noise
        rates = (X(:,transient+1:end) + 100)/20;
        if ss == 1
            rates = rates + 1*randn(size(rates))/sqrt(N);
        else
            rates = rates + .5*randn(1,size(rates,2))/sqrt(N);
        end

        % compute dimensionality
        Y = rates;
        Y = Y - mean(Y,2,'omitnan');
        Yz = sum(Y,2);
        if ~sum(Yz)
            continue;
        end
        indx = Yz ~= 0;
        Y = Y(indx,:);

        [~,~,latent] = pca(Y,'Centered','off');

        cve = cumsum(latent/sum(latent));
        II = find(cve > .8);
        if ~isempty(II)
            dim_r(ss) = II(1);
        end

        % compute stability
        pseudo_rates = nan(size(rates));
        for c = 1:N
            indx = 1:N;
            indx(c) = [];
            pseudo_rates(c,:) = mean(rates(indx,:),1);
        end

        windsz = 5;
        timewindow = [1 : windsz];
        X1 = pseudo_rates(:,timewindow);
        X2 = pseudo_rates(:,end-windsz+1 : end);

        mX1 = mean(X1,2);
        mX2 = mean(X2,2);
        stability = abs(mX1 - mX2) < mX1/100;
        no_change = stability == 1;

        pstable = 100*sum(no_change) / length(stability);
        stability_r(ss) = pstable;
        
    end
    dim_all(:,r) = dim_r;
    stability_all(:,r) = stability_r;
end

toc

%% plot
figure; 
bar(mean(stability_all,2)); hold on; ylabel('Percent stable networks');
errorbar(mean(stability_all,2), std(stability_all,[],2)/sqrt(Nr),'LineStyle','none','Color','k','LineWidth',2);
xticklabels({'Uncorrelated', 'Strongly correlated'}); box off;
xtickangle(45); set(gca,'fontsize',15);

figure; 
bar(mean(dim_all,2)); hold on; ylabel('Dimensionality');
errorbar(mean(dim_all,2), std(dim_all,[],2)/sqrt(Nr),'LineStyle','none','Color','k','LineWidth',2);
xticklabels({'Uncorrelated', 'Strongly correlated'}); box off;
xtickangle(35); set(gca,'fontsize',15);
