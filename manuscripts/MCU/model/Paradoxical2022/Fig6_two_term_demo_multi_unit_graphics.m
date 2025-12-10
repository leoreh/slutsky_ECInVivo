%called by Fig6_two_term_demo_multi_unit
%ssaray@ucla.edu
%dbuono@ucla.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOTTING %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(1)
     h1=figure('Renderer', 'painters', 'Position', [10 10 1200 1000]);
     set(gcf,'color','w');
     hold on
     
        subplot(5,2,1)
        plot(dt-EvokedOn*dt:dt:tmax*dt-EvokedOn*dt,hR(:,18),'color',[0 0.5 0],'linewidth',2,'ydatasource','hR(:,18)'); 
        line([dt-EvokedOn*dt tmax*dt-EvokedOn*dt],[ExSet ExSet],'color',[0 0.5 0],'linestyle',':','linewidth',2);
        ylim([0 20])
        xlim([dt-EvokedOn*dt 0.7])
        ylabel('E (Hz)','FontSize',16)
        set(gca,'FontSize',12)


        subplot(5,2,3)
        plot(dt-EvokedOn*dt:dt:tmax*dt-EvokedOn*dt,hR(:,40),'color',[0 0.5 0],'linewidth',2,'ydatasource','hR(:,40)');   
        line([dt-EvokedOn*dt tmax*dt-EvokedOn*dt],[ExSet ExSet],'color',[0 0.5 0],'linestyle',':','linewidth',2);
        ylim([0 20])
        xlim([dt-EvokedOn*dt 0.7])
        set(gca,'FontSize',12)

        
        subplot(5,2,5)
        plot(dt-EvokedOn*dt:dt:tmax*dt-EvokedOn*dt,hR(:,Ne+1),'color',[1 0 0],'linewidth',2,'ydatasource','hR(:,Ne+1)');  
        line([dt-EvokedOn*dt tmax*dt-EvokedOn*dt],[InhSet InhSet],'color',[1 0 0],'linestyle',':','linewidth',2);
        ylim([0 55])
        xlim([dt-EvokedOn*dt 0.7])
        ylabel('I (Hz)','FontSize',16)
        set(gca,'FontSize',12)


        subplot(5,2,7)
        plot(dt-EvokedOn*dt:dt:tmax*dt-EvokedOn*dt,hR(:,Ne+2),'color',[1 0 0],'linewidth',2,'ydatasource','hR(:,Ne+2)');
        line([dt-EvokedOn*dt tmax*dt-EvokedOn*dt],[InhSet InhSet],'color',[1 0 0],'linestyle',':','linewidth',2);
        xlabel('Time (sec)','Fontsize',16)
        ylim([0 55])
        xlim([dt-EvokedOn*dt 0.7])
        set(gca,'FontSize',12)

        
     subplot(5,2,[9 10])
     plot(NaN(nTrial,Ni),'color',[255/255,112/255,112/255],'ydatasource','trialhistFCaInh','linewidth',1); 
     hold on
     plot(NaN(nTrial,20),'color',[77/255 166/255 77/255],'ydatasource','trialhistFCa(:,1:20)','linewidth',1);
     line([1 nTrial],[ExSet ExSet],'color',[0 0.5 0],'linestyle','--','linewidth',2);
     line([1 nTrial],[InhSet InhSet],'color',[255/255,51/255,51/255],'linestyle','--','linewidth',2);
     ylabel('Mean E/I (Hz)','Fontsize',16)
     xlabel('Trials','Fontsize',16)
     ylim([0 35]) 
     set(gca,'FontSize',14)

        
        
     subplot(5,2,[2 4])
     plot(trialhistWEEp,trialhistWEIp,'-o','MarkerSize',3,'LineWidth',0.5,'color',[196/255,193/255,193/255],'xdatasource','trialhistWEEp','ydatasource','trialhistWEIp');
     hold on
     plot(WEEp,WEIp,'o','LineWidth',2,'color',[0 0.5 0],'xdatasource','WEEp','ydatasource','WEIp');
     xlabel('WEE')
     ylabel('WEI')
     set(gca,'FontSize',12)
     hold on
     x=WEEp;
     y1=WEIp;
     [R,p] = corr(x,y1,'rows','complete');
     P = polyfit(x,y1,1);
     yfit = P(1)*x+P(2);
     plot(x,yfit,'b','LineWidth',2,'xdatasource','x','ydatasource','yfit');
     
        
     subplot(5,2,[6 8])
     plot(trialhistWIEp,trialhistWIIp,'-o','MarkerSize',3,'LineWidth',0.5,'color',[196/255,193/255,193/255],'xdatasource','trialhistWIEp','ydatasource','trialhistWIIp');
     hold on
     plot(WIEp,WIIp,'o','color',[1 0 0],'LineWidth',2,'xdatasource','WIEp','ydatasource','WIIp');
     xlabel('WIE')
     ylabel('WII')
     set(gca,'FontSize',12)
     hold on
     ix=WIEp;
     iy1=WIIp;
     [iR,ip] = corr(ix,iy1,'rows','complete');
     iP = polyfit(ix,iy1,1);
     iyfit = iP(1)*ix+iP(2);
     plot(ix,iyfit,'b','LineWidth',2,'xdatasource','ix','ydatasource','iyfit');
end


