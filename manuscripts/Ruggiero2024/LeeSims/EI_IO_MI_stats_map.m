
%  This script simulates a balanced E/I network driven by a
%  piecewise-constant input signal, and computes the mutual information
%  between input and population-averaged firing-rate, and single-unit
%  firing-rates (Supplementary Figure 2).
%
%  The code loops over a range of input change-rates and a range of
%  networks, where the excitatory-to-excitatory connectivity strength is
%  varied. For each combination, several random realizations are simulated.
%
clear all; close all;
clc;
tic

plot_figs = false;

N = 1000;
Nr = 10;
NE = .8 * N;
NI = .2 * N;
K = N/10;

% Time params and variables
t_end = 5000;
dt = 1/10;
nsecs = t_end;
perttime = round((nsecs/2)/dt);
simtime = 0:dt:nsecs-dt;
simtime_len = length(simtime);
transient = simtime_len/4;
steady_state = simtime_len/2;
% Neural activation function
f = @(x,b) (x+b).*(x+b>0);

% Connectivity params
coupling = K/N;
J_IE = 10;
J_EI = -2.6;  J_II = -12;
J_exti = 1; J_exte = 1;
p = 1;

% Input and EE connectivity param ranges
freq_vals = [1,2,5,10,15,20,25,30,35,40,45,50,75,100];
Nfreq = length(freq_vals);
NJ_EE = 15;
J_EE_vals = linspace(.5,1.6,NJ_EE);

mfr_all = nan(NJ_EE,Nfreq,Nr);
mutual_all = nan(NJ_EE,Nfreq,Nr);
single_mutual_all = nan(NJ_EE,Nfreq,N,Nr);

for i=1:NJ_EE
    J_EE = J_EE_vals(i);

    for rr = 1:Nr
        
        % Build random connectivity matrix
        Aee = binornd(1,coupling,NE,NE);
        Aii = binornd(1,coupling,NI,NI);
        Aei = binornd(1,coupling,NE,NI);
        Aie = binornd(1,coupling,NI,NE);
        E_input_mask = binornd(1,p,NE,1);
        I_input_mask = binornd(1,p,NI,1);
        Jmat = [J_EE.*Aee J_EI.*Aei; J_IE.*Aie J_II.*Aii]./sqrt(K);

        x0 = rand(N,1);
        h0 = [J_exte.*ones(NE,1).*E_input_mask; J_exti.*ones(NI,1).*I_input_mask]*sqrt(K/p);
        b = 1;
        g = 1;

        % loop over input change rates
        for j=1:Nfreq

            freq = freq_vals(j);

            X = nan(N,simtime_len);
            R = nan(N,simtime_len);
            z_all = nan(simtime_len,1);

            X(:,1) = x0;
            x = x0;
            r = f(x,b);
            ti = 0;

            nhours = 1;
            windsz = round(nhours * (nsecs/dt)/20);

            n = [ones(NE,1)/N; ones(NI,1)/N] + 0*randn(N,1);
            m = [ones(NE,1); ones(NI,1)].*rand(N,1);

            bout_len = freq*10/dt;
            nbouts = floor(simtime_len / bout_len);
            for ii = 1:nbouts-1
                input(ii*bout_len + 1 : (ii+1)*bout_len ) = rand*ones(1,bout_len);
            end
            input(1:steady_state) = 1;

            % simulate network dynamics
            for ti = 1:simtime_len

                r = f(x,b);
                X(:,ti+1) = (1.0-dt)*x + g.*(Jmat*r + h0*input(ti))*dt;
                R(:,ti) = r;
                x = X(:,ti+1);

                z = n'*r;
                z_all(ti) = z;

            end
            %%%%%%%%%%%%%%%%%%%%

            R = R(:,transient+1:end);
            YE = R(1:NE,:);
            YI = R(NE+1:end,:);
            mean_FR = mean(R(:,1:steady_state-transient),2);

            z_ss = z_all(steady_state+1:end);
            input_ss = input(steady_state+1:end)';
            MFR = mean(mean_FR);
            if MFR > 1e5
                MFR = NaN;
            end
            mfr_all(i,j,rr) = MFR;
            mutual_all(i,j,rr) = mi(input_ss, z_ss);

            %%% Compute single-cell mutual information
            [~, fr_indx] = sort(mean_FR,'descend');
            for ni = 1:N
                single = R(fr_indx(ni),steady_state-transient+1:end)';
                single_mutual_all(i,j,ni,rr) = mi(single, z_ss);
            end
        end
    end
end
toc
%% plot

mmfr_all = mean(mfr_all,3,'omitnan');
mmutual_all = mean(mutual_all,3);
mm = mean(mmfr_all,2,'omitnan');
figure; imagesc(mm, 1./(10*freq_vals), mmutual_all');
set(gca,'YDir','normal','YScale','log');
ylabel('Input change rate'); xlabel('MFR'); colorbar;
title('Mutual Information'); 

%%
msingle_mutual_all = mean(single_mutual_all,4,'omitnan');
mmmfr_all = mean(mmfr_all,2);
m1 = squeeze(msingle_mutual_all(1,10,:));
m2 = squeeze(msingle_mutual_all(10,10,:));
figure;
plot(m1,'.','MarkerSize',25); hold on;
plot(m2,'.','MarkerSize',25); hold on;

%%

relative_mutual_pop = mmutual_all;
srelative_mutual_pop = std(relative_mutual_pop,[],2)/sqrt(size(relative_mutual_pop,2));
mmsingle_mutual_all = mean(msingle_mutual_all,3);
relative_mutual_single = mmsingle_mutual_all;
srelative_mutual_single = std(relative_mutual_single,[],2) / sqrt(size(relative_mutual_single,2));


figure;
[mmmfr_all,I] = sort(mmmfr_all,'ascend');
x = mmmfr_all';
y = mean(relative_mutual_pop(I),2)';
erry = srelative_mutual_pop(I)';
err = [y + erry; y - erry];
yyaxis left;
plot(x,y,'LineWidth',2,'Color','k'); hold on;
patch([x, fliplr(x)], [err(1,:), fliplr(err(2,:))], 'k', 'FaceAlpha',0.2, 'EdgeColor', 'none');
set(gca,'YColor','k');
y = mean(relative_mutual_single(I),2)';
erry = srelative_mutual_single(I)';
err = [y + erry; y - erry];
ylabel('MFR MI');
yyaxis right;
plot(x,y,'LineWidth',2,'Color',[0 .5 0]);
patch([x, fliplr(x)], [err(1,:), fliplr(err(2,:))],'g', 'FaceAlpha',0.2, 'EdgeColor', 'none')  
box off; axis square;
xlabel('MFR'); ylabel('Single-unit MI');
set(gca,'YColor','g');

relative_mutual_pop = mmutual_all;% ./ mmutual_all(:,end);
relative_mutual_single = mmsingle_mutual_all;% ./ mmsingle_mutual_all(:,end);
srelative_mutual_pop = std(relative_mutual_pop,[],1)/sqrt(size(relative_mutual_pop,1));
srelative_mutual_single = std(relative_mutual_single,[],1) / sqrt(size(relative_mutual_single,1));

figure;
x = 1./(10*freq_vals);
y = mean(relative_mutual_pop,1);
erry = srelative_mutual_pop;
err = [y + erry; y - erry];
yyaxis left;
plot(x,y,'LineWidth',2,'Color','k'); hold on;
patch([x, fliplr(x)], [err(1,:), fliplr(err(2,:))], 'k', 'FaceAlpha',0.2, 'EdgeColor', 'none');
set(gca,'YColor','k');
y = mean(relative_mutual_single,1);
erry = srelative_mutual_single;
err = [y + erry; y - erry];
ylabel('MFR MI');
yyaxis right;
plot(x,y,'LineWidth',2,'Color',[0 .5 0]);
patch([x, fliplr(x)], [err(1,:), fliplr(err(2,:))],'g', 'FaceAlpha',0.2, 'EdgeColor', 'none')  
box off; axis square;
set(gca,'xscale','log','YColor','g');
xlabel('Input change rate'); ylabel('Single-unit MI');
