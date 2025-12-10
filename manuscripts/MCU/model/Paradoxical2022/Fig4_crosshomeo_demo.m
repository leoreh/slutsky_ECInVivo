%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rate-based population model with cross-homeostatic ISN development rules 
% Based on Jercog et al (2017) https://elifesciences.org/articles/22425
% ssaray@ucla.edu
% dbuono@ucla.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SIMULATION SETTINGS %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

dt = 0.0001; %in seconds
tmax   = 2/dt; 
nTrial = 500;

rng(42) 
    
VIDEO =0;
HOMEOSTATIC_FLAG = 1;
GRAPHICS = 1;

learning_rule= 'cross_homeo'; 

%savetrials=[1,20,100,nTrial]; %NOTE: displayed trials numbers on the paper figure were rounded for simplicity. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NEURON PARAMETERS %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

F = @(x,gain,Thr) gain*max(0,x-Thr); %transfer function

thetaE = 4.8;
thetaI = 25;
gainE = 1;
gainI = 4;

Etau = 10/(dt*1000); %in miliseconds
Itau = 2/(dt*1000);  

Beta = 0;
tauA = 500;

E_MAX = 100; % Saturation of excitatory and inhibitory neurons
I_MAX = 250;

OUtau = 0.1; % Ornstein Uhlenbeck Noise
OUmu = 0;
OUsigma = 0.1; %sigma * sqrt(dt)

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

alpha = 0.0005; %learning rate 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INIT VARS %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

WEE = 2.1;   
WEI = 3;    
WIE = 4;    
WII = 1.5;   

WInit = [WEE WEI WIE WII];

hR = zeros(tmax,2); %history of Inh and Ex rate and adaptation

EvokedOn = 0.250/dt; %current injection to elicit a permanent Up
EvokedDur = 0.01/dt;
EvokedAmp = 7;

trialhistFCa      = NaN(nTrial,1);
trialhistFCaInh   = NaN(nTrial,1);
trialhistWEE      = NaN(nTrial,1);
trialhistWEI      = NaN(nTrial,1);
trialhistWIE      = NaN(nTrial,1);
trialhistWII      = NaN(nTrial,1);  
   
ExAvg = 0;
InhAvg = 0;
      
OUE = 0;
OUI = 0;

counter = 0;  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOTTING %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if GRAPHICS
    
    h = figure('Position',[1 1 1200 900]);
    set(gcf,'color','w');
    
    subplot(3,1,1)
    colormap = [ 255/255,51/255,51/255; 0 0.5 0; 0 0 1];
    set(gca,'colororder',colormap);
    hold on
    plot(dt-EvokedOn*dt:dt:tmax*dt-EvokedOn*dt,hR,'ydatasource','hR','linewidth',3);
    line([dt-EvokedOn*dt tmax*dt-EvokedOn*dt],[ExSet ExSet],'color',[0 0.5 0],'linestyle','--','linewidth',2);
    line([dt-EvokedOn*dt tmax*dt-EvokedOn*dt],[InhSet InhSet],'color',[255/255,51/255,51/255],'linestyle','--','linewidth',2);
    ylabel('E/I (Hz)')
    xlabel('Time (sec)')
    str = sprintf('Ex(green) Inh(red)');
    set(gca,'FontSize',20)
    set(gca,'linew',2)
    set(gca, 'box', 'off')
    ylim([0 20])
    xlim([dt-EvokedOn*dt tmax*dt-EvokedOn*dt])

    subplot(3,1,2)
    plot(zeros(nTrial,1),'color',[255/255,51/255,51/255],'ydatasource','trialhistFCaInh','linewidth',3);
    hold on
    plot(zeros(nTrial,1),'color',[0 0.5 0],'ydatasource','trialhistFCa','linewidth',3);
    line([1 nTrial],[ExSet ExSet],'color',[0 0.5 0],'linestyle','--','linewidth',2);
    line([1 nTrial],[InhSet InhSet],'color',[255/255,51/255,51/255],'linestyle','--','linewidth',2);
    ylabel('Mean E/I (Hz)')
    set(gca,'FontSize',20)
    set(gca,'linew',2)
    set(gca, 'box', 'off')
    ylim([0 20])

    subplot(3,1,3)
    plot(zeros(nTrial,1),'color',[18/255, 181/255, 143/255],'ydatasource','trialhistWEE','linewidth',3);
    hold on
    plot(zeros(nTrial,1),'color',[18/255, 181/255, 143/255],'ydatasource','trialhistWEI','linestyle',':','linewidth',3);
    plot(zeros(nTrial,1),'color',[217/255, 68/255, 220/255],'ydatasource','trialhistWIE','linewidth',3);
    plot(zeros(nTrial,1),'color',[217/255, 68/255, 220/255],'ydatasource','trialhistWII','linestyle',':','linewidth',3);
    xlim([0 nTrial])
    ylabel('Weights')
    xlabel('Trials')
    str = sprintf('g-=WEE g:=WEI m-=WIE m:=WII');
    legend('WEE','WEI','WIE','WII')
    set(gca,'FontSize',20)
    set(gca,'linew',2)
    set(gca, 'box', 'off')
    legend('WEE','WEI','WIE','WII','LineWidth',1)
    ylim([0 6.5])


end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SIMULATION %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic

    for trial=1:nTrial
      
      
      E = 0; %init rates 
      I = 0;
      a = 0;
      
      
      fCa = zeros(1,tmax);  %instantaneous fast Ca sensor
      fCaInh = zeros(1,tmax);  %instantaneous fast Ca sensor
      
      evoked = zeros(1,tmax);
      evoked(EvokedOn:EvokedOn+EvokedDur)=EvokedAmp;

      
      hR = zeros(tmax,2); %history of Inh and Ex rate and adaptation
      
      for t=1:tmax
         
         
         
         OUE = OUE + OUtau*(OUmu-OUE) + OUsigma*randn; %Ornstein-Uhlenbeck Noise for excitatory unit
         OUI = OUI + OUtau*(OUmu-OUI) + OUsigma*randn; %Ornstein-Uhlenbeck Noise for inhibitory unit
         
         E = E + (-E + F(WEE*E - WEI*I - a + evoked(t) + OUE ,gainE,thetaE) )/Etau;
         I = I + (-I + F(WIE*E - WII*I + OUI,gainI,thetaI) )/Itau;
          
         if E>E_MAX; E = E_MAX; end
         if I>I_MAX; I = I_MAX; end
          
         %Ca Sensors
         fCa(:,t) = E;
         fCaInh(:,t) = I;
         
         hR(t,:) = [I E];
         
      end
      
      
      %% HOMEOSTASIS
      
         ExAvg    = ExAvg  + (-ExAvg  + mean(fCa((end-0.5/dt):end)))/tau_trial; %integrate at the end of trial to avoid evoked
         InhAvg   = InhAvg + (-InhAvg + mean(fCaInh((end-0.5/dt):end)))/tau_trial;
      
         
      if HOMEOSTATIC_FLAG
          
          EAvg =  max(1,ExAvg); %Average activity is rectified, for trials that start with 0 rate (development settings), otherwise weights would never move 
          IAvg = max(1,InhAvg); 
          
            newWEE = WEE + alpha*EAvg*(InhSet-IAvg);
            newWEI = WEI - alpha*IAvg*(InhSet-IAvg);
            newWIE = WIE - alpha*EAvg*(ExSet-EAvg); 
            newWII = WII + alpha*IAvg*(ExSet-EAvg);
            
         WEE = newWEE; WEI = newWEI; WIE = newWIE; WII = newWII;
         
         if WEE<WEE_MIN; WEE = WEE_MIN; end
         if WEI<WEI_MIN; WEI = WEI_MIN; end
         if WIE<WIE_MIN; WIE = WIE_MIN; end
         if WII<WII_MIN; WII = WII_MIN; end
      
      end
      
      
      trialhistFCa(trial) = ExAvg;
      trialhistFCaInh(trial) = InhAvg;
      trialhistWEE(trial) = WEE;
      trialhistWEI(trial) = WEI;
      trialhistWIE(trial) = WIE;
      trialhistWII(trial) = WII;
      
      if GRAPHICS %&& rem(trial,10)==0
      refreshdata 
      drawnow
%         if ismember(trial,savetrials)
%             saveas(gcf,[learning_rule,'_demo','/temp','Ex',num2str(ExSet),'Inh',num2str(InhSet),'Trial',num2str(trial)],'tiff') 
%             saveas(gcf,[learning_rule,'_demo','/temp','Ex',num2str(ExSet),'Inh',num2str(InhSet),'Trial',num2str(trial)],'epsc') 
%         end        
      end
      if VIDEO
      counter = counter+1;
      frames(counter) = getframe(h); 
      end
    end
  
   if VIDEO
    video = VideoWriter([learning_rule,'_demo','/exp',num2str(exp)]);
    video.FrameRate = 10; 
    open(video)
    writeVideo(video,frames);
    close(video)
    clear video
    clear frames
   end
toc
%save([learning_rule,'_demo','/temp','Ex',num2str(ExSet),'Inh',num2str(InhSet)])
