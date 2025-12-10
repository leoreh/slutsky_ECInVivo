%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ONLINE Rate-based population model with cross-homeostatic ISN development rules 
% Based on Jercog et al (2017) https://elifesciences.org/articles/22425
% ssaray@ucla.edu
% dbuono@ucla.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SIMULATION SETTINGS %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

rng(42) 

dt = 0.0001; %IN SECONDS
tmax   = 100/dt; 


GRAPHICS = 1;
HOMEOSTATIC_FLAG = 1;
VIDEO = 0; 

learning_rule = 'cross-homeo_ONLINE';

%% NETWORK PARAMETERS

F = @(x,gain,Thr) gain*max(0,x-Thr);

WEE = 2.1; 
WEI = 3;    
WIE = 4;  
WII = 2;  

thetaE = 4.8;
thetaI = 25;
gainE = 1;
gainI = 4;

Etau = 10/(dt*1000); %10
Itau = 2/(dt*1000);  %2

Beta = 0.0;
tauA = 500;

E_MAX = 100; % Saturation of excitatory and inhibitory neurons
I_MAX = 250; 

WEI_MIN = 0.1;
WEE_MIN = 0.1;
WII_MIN = 0.1;
WIE_MIN = 0.1;

% Ornstein Uhlenbeck Noise
OUtau = 1/10;
OUmu = 0;
OUsigma = 0.1; %sigma * sqrt(dt)

%% Evoked: current injection to elicit a permanent Up
% drawn from a poisson distribution 

poisson_rate = 0.1; %in Hz
poisson_dur = tmax*dt; % duration of the trial in seconds
number_trains = 1; % number of trains to generate 

[spikevec, tVec] = poissonSpikeGen(poisson_rate, poisson_dur, number_trains, dt );

EvokedOn = double(spikevec);%
EvokedDur = 0.01/dt;%
EvokedAmp = 7; %

spkidx = find(EvokedOn);
spkidx = [0.2/dt,spkidx]; %we force one evoked event at the beginning of the simulation
spkdt=spkidx+EvokedDur;
      
evoked = zeros(1,tmax);
for m=1:length(spkidx)
evoked(spkidx(m):spkdt(m))=EvokedAmp;
end

% figure
% plot(evoked)


%% PLASTICITY STUFF

ExSet = 5; 
InhSet = 14; 

alpha = 0.0000005; %learning rate 0.0005

tauFCa=1000/(dt*1000); %time constant of calcium sensors in ms
tauFCaInh=1000/(dt*1000); 


%% Init Variables

WInit = [WEE WEI WIE WII];
   
ExAvg = 0;
InhAvg = 0;
      
counter = 0;        
      
E = 0;
I = 0;
a = 0;
      
FCa =0; %Ca sensor for E 
FCaInh =0; %Ca sensor for I

fCa = NaN(1,tmax);  %Ca sensor time vec
fCaInh = NaN(1,tmax); 

histWEE = NaN(1,tmax);
histWEI = NaN(1,tmax);
histWIE = NaN(1,tmax);
histWII = NaN(1,tmax);

hR = NaN(tmax,2); %history of Inh and Ex rate 

OUE = 0;
OUI = 0;

%% GRAPHICS

if GRAPHICS

h = figure('Position',[1 1 600 800]);
        set(gcf,'color','w');

subplot(3,1,1)
colormap = [ 255/255,51/255,51/255; 0 0.5 0; 0 0 1];
set(gca,'colororder',colormap);
hold on
plot(dt:dt:tmax*dt,hR,'ydatasource','hR','linewidth',3);
line([dt tmax*dt],[ExSet ExSet],'color',[0 0.5 0],'linestyle','--','linewidth',2);
line([dt tmax*dt],[InhSet InhSet],'color',[255/255,51/255,51/255],'linestyle','--','linewidth',2);
ylabel('E/I (Hz)')
xlabel('Time (sec)')
str = sprintf('Ex(green) Inh(red)');
set(gca,'FontSize',20)
set(gca,'linew',2)
set(gca, 'box', 'off')
ylim([0 20])

subplot(3,1,2)
hold on
plot(dt:dt:tmax*dt,fCa,'ydatasource','fCa','linewidth',3,'color',[0 0.5 0]);
plot(dt:dt:tmax*dt,fCaInh,'ydatasource','fCaInh','linewidth',3,'color',[255/255,51/255,51/255]);
line([dt tmax*dt],[ExSet ExSet],'color',[0 0.5 0],'linestyle','--','linewidth',2);
line([dt tmax*dt],[InhSet InhSet],'color',[255/255,51/255,51/255],'linestyle','--','linewidth',2);
ylabel('Calcium E/I (Hz)')
xlabel('Time (sec)')
str = sprintf('Ex(green) Inh(red)');
set(gca,'FontSize',20)
set(gca,'linew',2)
set(gca, 'box', 'off')
ylim([0 20])

subplot(3,1,3)
plot(dt:dt:tmax*dt,histWEE,'color',[18/255, 181/255, 143/255],'ydatasource','histWEE','linewidth',3);
hold on
plot(dt:dt:tmax*dt,histWEI,'color',[18/255, 181/255, 143/255],'ydatasource','histWEI','linestyle',':','linewidth',3);
plot(dt:dt:tmax*dt,histWIE,'color',[217/255, 68/255, 220/255],'ydatasource','histWIE','linewidth',3);
plot(dt:dt:tmax*dt,histWII,'color',[217/255, 68/255, 220/255],'ydatasource','histWII','linestyle',':','linewidth',3);
line([dt tmax*dt],[1 1],'color','w','linestyle','--','linewidth',2);
ylabel('Weights')
xlabel('Time (sec)')
legend('W_{EE}','W_{EI}','W_{IE}','W_{II}','LineWidth',1)
set(gca,'FontSize',20)
set(gca,'linew',2)
set(gca, 'box', 'off')
ylim([0 6]) 


end

%% SIMULATION 

tic                  
      for t=1:tmax
         
         
         OUE = OUE + OUtau*(OUmu-OUE) + OUsigma*randn; %Ornstein-Uhlenbeck Noise for excitatory unit
         OUI = OUI + OUtau*(OUmu-OUI) + OUsigma*randn; %Ornstein-Uhlenbeck Noise for inhibitory unit
         
         E = E + (-E + F(WEE*E - WEI*I - a + evoked(t) + OUE ,gainE,thetaE) )/Etau;
         I = I + (-I + F(WIE*E - WII*I + OUI,gainI,thetaI) )/Itau;
          
         if E>E_MAX; E = E_MAX; end
         if I>I_MAX; I = I_MAX; end
          
         % Ca Sensors
         FCa = FCa + (-FCa + E)/tauFCa; 
         FCaInh = FCaInh + (-FCaInh + I)/tauFCaInh; 

         
         fCa(:,t) = FCa;
         fCaInh(:,t) = FCaInh;
         
         hR(t,:) = [I E];
         
         histWEE(:,t) = WEE;
         histWEI(:,t) = WEI;
         histWIE(:,t) = WIE;
         histWII(:,t) = WII;
      
         
        %% HOMEOSTASIS
      
         
          if HOMEOSTATIC_FLAG 

              EAvg =  FCa; 
              IAvg = FCaInh; 

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
          
             if GRAPHICS && rem(t,1/dt)==0
              refreshdata %refresh graphics
              drawnow
            end
         
      end
        
      
   if VIDEO
    counter = counter+1;
    frames(counter) = getframe(h); 
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

% saveas(gcf,'historyEIW','jpg')
% saveas(gcf,'historyEIW','svg')
    

%% PLOTTING

figure('Position',[1 1 600 800]);
set(gcf,'color','w');

subplot(3,1,1)
colormap = [ 255/255,51/255,51/255; 0 0.5 0; 0 0 1];
set(gca,'colororder',colormap);
hold on
plot(dt:dt:tmax*dt,hR,'ydatasource','hR','linewidth',3);
line([dt tmax*dt],[ExSet ExSet],'color',[0 0.5 0],'linestyle','--','linewidth',2);
line([dt tmax*dt],[InhSet InhSet],'color',[255/255,51/255,51/255],'linestyle','--','linewidth',2);
ylabel('E/I (Hz)')
xlabel('Time (sec)')
str = sprintf('Ex(green) Inh(red)');
set(gca,'FontSize',20)
set(gca,'linew',2)
set(gca, 'box', 'off')
ylim([0 20])
xlim([0 1])

subplot(3,1,2)
colormap = [ 255/255,51/255,51/255; 0 0.5 0; 0 0 1];
set(gca,'colororder',colormap);
hold on
plot(dt:dt:tmax*dt,hR,'ydatasource','hR','linewidth',3);
line([dt tmax*dt],[ExSet ExSet],'color',[0 0.5 0],'linestyle','--','linewidth',2);
line([dt tmax*dt],[InhSet InhSet],'color',[255/255,51/255,51/255],'linestyle','--','linewidth',2);
ylabel('E/I (Hz)')
xlabel('Time (sec)')
str = sprintf('Ex(green) Inh(red)');
set(gca,'FontSize',20)
set(gca,'linew',2)
set(gca, 'box', 'off')
ylim([0 20])
xlim([20 30])

subplot(3,1,3)
colormap = [ 255/255,51/255,51/255; 0 0.5 0; 0 0 1];
set(gca,'colororder',colormap);
hold on
plot(dt:dt:tmax*dt,hR,'ydatasource','hR','linewidth',3);
line([dt tmax*dt],[ExSet ExSet],'color',[0 0.5 0],'linestyle','--','linewidth',2);
line([dt tmax*dt],[InhSet InhSet],'color',[255/255,51/255,51/255],'linestyle','--','linewidth',2);
ylabel('E/I (Hz)')
xlabel('Time (sec)')
str = sprintf('Ex(green) Inh(red)');
set(gca,'FontSize',20)
set(gca,'linew',2)
set(gca, 'box', 'off')
ylim([0 20])
xlim([99 100])
% saveas(gcf,'historyEIWzoom','jpg')
% saveas(gcf,'historyEIWzoom','svg')
