%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multi-unit rate-based model with cross-homeostatic ISN development rules
% Based on Jercog et al (2017) https://elifesciences.org/articles/22425
% ssaray@ucla.edu
% dbuono@ucla.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SIMULATION SETTINGS %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clear all
% close all

dt = 0.0001; %IN SECONDS
tmax   = 2/dt; 
nTrial = 250; 

rng(44)
    
GRAPHICS = 1;
VIDEO =0;
HOMEOSTATIC_FLAG = 1;


learning_rule= 'cross_homeo'; 

%savetrials = [1,2,5,40,nTrial]; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NEURON PARAMETERS %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F = @(x,gain,Thr) gain*max(0,x-Thr);

Ne = 80;
Ni = 20;

thetaE = 4.8;
thetaI = 25;
gainE = 1;
gainI = 4;

Etau = 10/(dt*1000); 
Itau = 2/(dt*1000);  

Beta = 0;
tauA = 500;

E_MAX = 100;
I_MAX = 250;

OUtau = 0.1; %Ornstein-Uhlenbeck Noise 
OUmu = 0;
OUsigma = 0.1; %sigma * sqrt(dt)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INIT WEIGHT MATRIX %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sigmaw = 0.04; %normally distributed weights with meanw and sigmaw
meanw = 0.1; 

WEE = meanw +randn(Ne,Ne)*sigmaw; 
WEE = WEE - diag(diag(WEE));

WEI = meanw +randn(Ne,Ni)*sigmaw;

WIE= meanw +randn(Ni,Ne)*sigmaw;

WII = meanw +randn(Ni,Ni)*sigmaw;
WII = WII - diag(diag(WII));


%%
Winit = [WEE,WEI;WIE,WII];

WEEp = sum(WEE,2); %sum of presynaptic weights
WEIp = sum(WEI,2);
WIEp = sum(WIE,2);
WIIp = sum(WII,2);

initWEEp = WEEp;
initWEIp = WEIp;
initWIEp = WIEp;
initWIIp = WIIp;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLASTICITY %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ExSet = 5; %setpoint for excitatory neurons
InhSet = 14; %setpoint for inhibitory neurons

tau_trial = 2 ;

WEI_MIN = 0.1;
WEE_MIN = 0.1;
WII_MIN = 0.1;
WIE_MIN = 0.1;

alpha = 0.00002; %learning rate 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INIT VARS %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
trialhistFCa      = NaN(nTrial,Ne);
trialhistFCaInh   = NaN(nTrial,Ni);
trialhistWEE      = NaN(nTrial,Ne);
trialhistWEI      = NaN(nTrial,Ne);
trialhistWIE      = NaN(nTrial,Ni);
trialhistWII      = NaN(nTrial,Ni);  
trialhistWEEp      = NaN(nTrial,Ne);
trialhistWEIp      = NaN(nTrial,Ne);
trialhistWIEp      = NaN(nTrial,Ni);
trialhistWIIp      = NaN(nTrial,Ni);  

ExAvg = zeros(Ne,1);
InhAvg = zeros(Ni,1);

hR = zeros(tmax,Ne+Ni); %history of Inh and Ex rate 


EvokedOn = 0.250/dt; %Evoked current
EvokedDur = 0.01/dt; 
EvokedAmp = 7; 

counter = 0;  

OUE = zeros(Ne,1);
OUI = zeros(Ni,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOTTING %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if GRAPHICS
    Fig5_crosshomeo_demo_multi_unit_graphics
end

if VIDEO
      counter = counter+1;
      frames(counter) = getframe(h1); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SIMULATION %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for trial=1:nTrial    
      
      fCa = zeros(Ne,tmax);  %instantaneous fast Ca sensor, integrates the firing rate of E
      fCaInh = zeros(Ni,tmax);  %instantaneous fast Ca sensor, integrates the firing rate of I
      
      hR = zeros(tmax,Ne+Ni); %history of Inh and Ex rate 
      
     E = zeros(Ne,1);
     I = zeros(Ni,1);
     a = zeros(Ne,1);
      
           
      evoked = zeros(1,tmax);
      evoked(EvokedOn:EvokedOn+EvokedDur)=EvokedAmp;
      
%%
      for t=1:tmax        
         
         OUE = OUE + OUtau*(OUmu-OUE) + OUsigma*randn(Ne,1); %Ornstein-Uhlenbeck Noise for excitatory unit
         OUI = OUI + OUtau*(OUmu-OUI) + OUsigma*randn(Ni,1); %Ornstein-Uhlenbeck Noise for inhibitory unit
         
         E = E + (-E + F(WEE*E - WEI*I - a + evoked(t) + OUE,gainE,thetaE) )/Etau;
         I = I + (-I + F(WIE*E - WII*I + OUI,gainI,thetaI) )/Itau;
         
         a = a + (-a + Beta*E)/tauA; %ADAPTATION 
         
         Emaxvec = E>E_MAX; E(Emaxvec) = E_MAX; %Neurons have a saturation of their rates
         Imaxvec = I>I_MAX; I(Imaxvec) = I_MAX;

         
         hR(t,:) = [E' I'];
         
         % Ex Ca Sensors
         fCa(:,t) = E;
         fCaInh(:,t) = I;       
         
      end
      
         ExAvg    = ExAvg  + (-ExAvg  + mean(fCa(:,end-0.5/dt:end),2))/tau_trial; %we average at the end of the trial to avoid evoked
         InhAvg   = InhAvg + (-InhAvg + mean(fCaInh(:,end-0.5/dt:end),2))/tau_trial;
         
      
    %%  
    
      WEEp = sum(WEE,2); %update sum of presynaptic weights for plot
      WEIp = sum(WEI,2);
      WIEp = sum(WIE,2);
      WIIp = sum(WII,2);
     
      trialhistFCa(trial,:) = ExAvg;
      trialhistFCaInh(trial,:) = InhAvg;
      trialhistWEE(trial,:) = WEE(:,end);      
      trialhistWEI(trial,:) = WEI(:,end);     
      trialhistWIE(trial,:) = WIE(:,end);   
      trialhistWII(trial,:) = WII(:,end);            
      trialhistWEEp(trial,:) = WEEp;
      trialhistWEIp(trial,:) = WEIp;
      trialhistWIEp(trial,:) = WIEp;
      trialhistWIIp(trial,:) = WIIp;
      
      x=WEEp; %find correlation between E - I weights
      y1=WEIp;
      [R,p] = corr(x,y1,'rows','complete');
      P = polyfit(x,y1,1);
      yfit = P(1)*x+P(2);
      ix=WIEp;
      iy1=WIIp;
      [iR,ip] = corr(ix,iy1,'rows','complete');
      iP = polyfit(ix,iy1,1);
      iyfit = iP(1)*ix+iP(2);
      
      
      if HOMEOSTATIC_FLAG  

            EAvg =  max(1,ExAvg); %Average activity is rectified, for trials that start with 0 rate (development settings), otherwise weights would never move
            IAvg = max(1,InhAvg);
            
            
            newWEE = WEE + alpha*EAvg'*sum(InhSet-IAvg)/Ni;
            newWEE = newWEE - diag(diag(newWEE));
            
            newWEI = WEI - alpha*IAvg'*sum(InhSet-IAvg)/Ni;

            newWIE = WIE - alpha*EAvg'*sum(ExSet-EAvg)/Ne; 
         
            newWII = WII + alpha*IAvg'*sum(ExSet-EAvg)/Ne; 
            newWII = newWII - diag(diag(newWII));
                   
         
        WEE = newWEE; WEI = newWEI; WIE = newWIE; WII = newWII;
         
         %If weights fall below a minimum or are NaN set to minimum
         
         WEE(WEE<WEE_MIN/(Ne-1))=WEE_MIN/(Ne-1);
         WEE(isnan(WEE))=WEE_MIN/(Ne-1);
         WEI(WEI<WEI_MIN/Ni)=WEI_MIN/Ni;
         WEI(isnan(WEI))=WEI_MIN/Ni;
         WIE(WIE<WIE_MIN/Ne)=WIE_MIN/Ne;
         WIE(isnan(WIE))=WIE_MIN/Ne;
         WII(WII<WII_MIN/(Ni-1))=WII_MIN/(Ni-1);
         WII(isnan(WII))=WII_MIN/(Ni-1);
         WEE = WEE - diag(diag(WEE));
         WII = WII - diag(diag(WII));
            

      end
     %%      
      
      if GRAPHICS %&& rem(trial,10)==0
      refreshdata(h1) 
      drawnow   
      
%           if ismember(trial,savetrials)
%               saveas(gcf,['trial',num2str(trial)],'jpg')
%               saveas(gcf,['trial',num2str(trial)],'svg')
% 
%           end
      
      end
      
      if VIDEO
      counter = counter+1;
      frames(counter) = getframe(h1); 
      end
      
      
 end
   
 Wend = [WEE,WEI;WIE,WII];  
 %save(['METApopmodel_',learning_rule,'Ex',num2str(ExSet),'Inh',num2str(InhSet),'.mat'])
 
   if VIDEO
    video = VideoWriter([learning_rule,'_MultiUnit_demo']);
    video.FrameRate = 2;
    open(video)
    writeVideo(video,frames);
    close(video)
    clear video
    clear frames
   end
   
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% final PLOTTING %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if GRAPHICS   

%% Weight matrix winit vs wend

Winit2 = Winit .* [ones(Ne+Ni,Ne),ones(Ne+Ni,Ni)*-1] ;     % make inhibitory weights negative

Wend2 = Wend .* [ones(Ne+Ni,Ne),ones(Ne+Ni,Ni)*-1] ;

figure('Position', [10 10 1300 500]);     
     subplot(1,2,1)
     imagesc(Winit2)    
     cMap = getPyPlot_cMap('coolwarm',128);
     colormap(cMap);
     colorbar
     caxis([-0.25 0.25])
     xlabel('pre')
     ylabel('post')
     title('pre training')    
     subplot(1,2,2)
    imagesc(Wend2)
    colorbar
    caxis([-0.25 0.25])
    xlabel('pre')
    ylabel('post')
    title('post training')
%     saveas(gcf,'Winitendcolor','jpg')
%     saveas(gcf,'Winitendcolor','svg')

    %% Histogram winit vs wend
    figure
    subplot(2,2,1)
    histogram(Winit(1:Ne,1:Ne))
    hold on
    histogram(WEE)
    xlabel('WEE')
    subplot(2,2,2)
    histogram(Winit(1:Ne,Ne+1:end))
    hold on
    histogram(WEI)
    legend('pre','post')
    xlabel('WEI')
    subplot(2,2,3)
    histogram(Winit(Ne+1:end,1:Ne))
    hold on
    histogram(WIE)
    xlabel('WIE')
    subplot(2,2,4)
    histogram(Winit(Ne+1:end,Ne+1:end))
    hold on
    histogram(WII)
    xlabel('WII')
%     saveas(gcf,'Winitendhistogram','jpg')
    
    %% E-I init
     
    figure('Position', [10 10 1200 500]);     
     subplot(1,2,1)
     scatter(initWEEp,initWEIp,[60],'filled','MarkerFaceColor',[0 0.5 0],'MarkerFaceAlpha',0.6)
     xlabel('WEE init')
     ylabel('WEI init')
     set(gca,'FontSize',20)
    set(findobj(gca,'type','line'),'linew',3)
    set(gca,'linew',4)
    set(gca, 'box', 'off')
    hold on
    x=initWEEp;
    y1=initWEIp;
    [R,p] = corr(x,y1,'rows','complete');
    P = polyfit(x,y1,1);
    yfit = P(1)*x+P(2);
    plot(x,yfit,'b','LineWidth',2,'color',[0, 0.4470, 0.7410]);
    title(['R = ',num2str(R),'p = ',num2str(p)],'FontSize',15)

    subplot(1,2,2)
     scatter(initWIEp,initWIIp,[60],'filled','MarkerFaceColor',[1 0 0],'MarkerFaceAlpha',0.6)
     xlabel('WIE init')
     ylabel('WII init')
     set(gca,'FontSize',20)
      xlim([7.4 8.6])
      ylim([1.4 2.3])
    set(findobj(gca,'type','line'),'linew',3)
    set(gca,'linew',4)
    set(gca, 'box', 'off')
    hold on
    x=initWIEp;
    y1=initWIIp;
    [R,p] = corr(x,y1,'rows','complete');
    P = polyfit(x,y1,1);
    yfit = P(1)*x+P(2);
    plot(x,yfit,'b','LineWidth',2,'color',[0, 0.4470, 0.7410]);
    title(['R = ',num2str(R),'p = ',num2str(p)],'FontSize',15)
%     saveas(gcf,'preWinit','jpg')
%     saveas(gcf,'preWinit','svg')


  
    %% E-I end
    figure('Position', [10 10 1200 500]);     
     subplot(1,2,1)
     scatter(WEEp,WEIp,[60],'filled','MarkerFaceColor',[0 0.5 0],'MarkerFaceAlpha',0.6)
     xlabel('WEE final')
     ylabel('WEI final')
     set(gca,'FontSize',20)
    set(findobj(gca,'type','line'),'linew',3)
    set(gca,'linew',4)
    set(gca, 'box', 'off')
    hold on
    x=WEEp;
    y1=WEIp;
    [R,p] = corr(x,y1,'rows','complete');
    P = polyfit(x,y1,1);
    yfit = P(1)*x+P(2);
    plot(x,yfit,'b','LineWidth',2,'color',[0, 0.4470, 0.7410]);
    title(['R = ',num2str(R),'p = ',num2str(p)],'FontSize',15)

    subplot(1,2,2)
    scatter(WIEp,WIIp,[60],'filled','MarkerFaceColor',[1 0 0],'MarkerFaceAlpha',0.6)
     xlabel('WIE final')
     ylabel('WII final')
     set(gca,'FontSize',20)
    set(findobj(gca,'type','line'),'linew',3)
    set(gca,'linew',4)
    set(gca, 'box', 'off')
     xlim([7.8 9.2])
     ylim([0.9 1.8])
    hold on
    x=WIEp;
    y1=WIIp;
    [R,p] = corr(x,y1,'rows','complete');
    P = polyfit(x,y1,1);
    yfit = P(1)*x+P(2);
    plot(x,yfit,'b','LineWidth',2,'color',[0, 0.4470, 0.7410]);
    title(['R = ',num2str(R),'p = ',num2str(p)],'FontSize',15)
%     saveas(gcf,'preW','jpg')
%     saveas(gcf,'preW','svg')

   
%       %% Currents end
%     figure('Position', [10 10 1200 500]);     
%      subplot(1,2,1)
%      plot(WEE*E,WEI*I,'o','color',[0 0.5 0])
%      xlabel('WEE*E')
%      ylabel('WEI*I')
%      set(gca,'FontSize',20)
%     set(findobj(gca,'type','line'),'linew',3)
%     set(gca,'linew',4)
%     set(gca, 'box', 'off')
%     hold on
%     x=WEE*E;
%     y1=WEI*I;
%     [R,p] = corr(x,y1,'rows','complete');
%     P = polyfit(x,y1,1);
%     yfit = P(1)*x+P(2);
%     plot(x,yfit,'b','LineWidth',2);
%     title(['R = ',num2str(R),'p = ',num2str(p)],'FontSize',15)
% 
%     subplot(1,2,2)
%      plot(WIE*E,WII*I,'o','color',[1 0 0])
%      xlabel('WIE*E')
%      ylabel('WII*I')
%      set(gca,'FontSize',20)
%     set(findobj(gca,'type','line'),'linew',3)
%     set(gca,'linew',4)
%     set(gca, 'box', 'off')
%     hold on
%     x=WIE*E;
%     y1=WII*I;
%     [R,p] = corr(x,y1,'rows','complete');
%     P = polyfit(x,y1,1);
%     yfit = P(1)*x+P(2);
%     plot(x,yfit,'b','LineWidth',2);
%     title(['R = ',num2str(R),'p = ',num2str(p)],'FontSize',15)
%     saveas(gcf,'currents','jpg')
%     
%     
%     %% Presynaptic vs postsynaptic weights 
%     
%       WEEpo = sum(WEE,1); %update sum of postsynaptic weights 
%       WEIpo = sum(WEI,1);
%       WIEpo = sum(WIE,1);
%       WIIpo = sum(WII,1);
%       
%     figure('Position', [10 10 1200 500]);     
%      subplot(1,2,1)
%      plot(WEEp,WEEpo,'o','color',[0 0.5 0])
%      xlabel('WEEp')
%      ylabel('WEEpo')
%      set(gca,'FontSize',20)
%     set(findobj(gca,'type','line'),'linew',3)
%     set(gca,'linew',4)
%     set(gca, 'box', 'off')
%     hold on
%     x=WEEp;
%     y1=WEEpo';
%     [R,p] = corr(x,y1,'rows','complete');
%     P = polyfit(x,y1,1);
%     yfit = P(1)*x+P(2);
%     plot(x,yfit,'b','LineWidth',2);
%     title(['R = ',num2str(R),'p = ',num2str(p)],'FontSize',15)
% 
%     subplot(1,2,2)
%      plot(WIEp,WEIpo,'o','color',[1 0 0])
%      xlabel('WIEp')
%      ylabel('WEIpo')
%      set(gca,'FontSize',20)
%     set(findobj(gca,'type','line'),'linew',3)
%     set(gca,'linew',4)
%     set(gca, 'box', 'off')
%     hold on
%     x=WIEp;
%     y1=WEIpo';
%     [R,p] = corr(x,y1,'rows','complete');
%     P = polyfit(x,y1,1);
%     yfit = P(1)*x+P(2);
%     plot(x,yfit,'b','LineWidth',2);
%     title(['R = ',num2str(R),'p = ',num2str(p)],'FontSize',15)
%     saveas(gcf,'outputW','jpg')


end

 