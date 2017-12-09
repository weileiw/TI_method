clear all
close all
load M3d.mat
load grid.mat
p.bf = 0.95;
figure(1)
dpa = 365; 
n234 = 10.5/dpa; n230 = 9.19e-6/dpa;
U8   = 2400*ones(24,1);
U4   = 2760*ones(24,1);
d0  = 0.001*ones(24,1);
d4  = 2.4*ones(24,1);
TH4 = 0*ones(24,1); 	
TH0 = 0*ones(24,1);
Th4  = 0*ones(24,1);	
Th0  = 0*ones(24,1);

TH4_b = TH4;
TH0_b = TH0;
TH4_b(1) = 2.5;
TH0_b(1) = 0.001;
k1 = 3/dpa;	  % aggregation
k2 = 150/dpa; % disagregation	
a1 = 0.5/dpa; % adsorption
a2 = 2/dpa;   % desorption
d1 = 1/dpa;   % remineralization

t = 0;	
dt = 1/96;
nstep = 10*365*24/dt;
for ik = 1:nstep

    dd4dt = U8*n234-d4*(a1+n234)+Th4*(a2+d1);
    dTh4dt = -Th4*(k1+d1+a2+n234)+d4*a1+TH4*k2;
    dTH4dt = -TH4*(k2+n234)+Th4*k1-100*(TH4+TH4_b-TH4([2:24,24])./grd.dzt);
    
    dd0dt = U4*n230-d0*(a1+n230)+Th0*(a2+d1); 
    dTh0dt = -Th0*(k1+d1+a2+n230)+d0*a1+TH0*k2;
    dTH0dt = -TH0*(k2+n230)+Th0*k1-100*(TH0+TH0_b-TH0([2:24,24])./grd.dzt);
    
    d4 = d4+dd4dt*dt;
    d0 = d0+dd0dt*dt;
    Th4 = Th4+dTh4dt*dt;
    Th0 = Th0+dTh0dt*dt;
    TH4  = TH4+dTH4dt*dt;
    TH0  = TH0+dTH0dt*dt;
    t = t+dt;
    
    if mod(ik,94*4) == 0
        
        set(0,'CurrentFigure',1);
        subplot(3,2,1)
        plot(t,d4(24),'rs'); hold on; drawnow
        ylabel('[^2^3^4Th_d] (dpm/m^3)')
        xlabel('time, (day)')
        
        subplot(3,2,2)
        plot(t,d0(24),'rs');hold on; drawnow
        ylabel('[^2^3^0Th_d (dpm/m^3)')
        xlabel('time, (day)')
        
        subplot(3,2,3)
        plot(t,TH4(24),'rs');hold on; drawnow
        ylabel('[^2^3^4Th_f] (dpm/m^3)')
        xlabel('time, (day)')
        
        subplot(3,2,4)
        plot(t,TH0(24),'rs'); hold on;drawnow
        ylabel('[^2^3^0Th_f] (dpm/m^3)')
        xlabel('time, (day)')
        
        subplot(3,2,5)
        plot(t,Th4(24),'rs');hold on;
        ylabel('[^2^3^4Th_s] (dpm/m^3)')
        xlabel('time, (day)')
        
        subplot(3,2,6)
        plot(t,Th0(24),'rs');hold on;drawnow
        ylabel('[^2^3^0Th_s] (dpm/m^3)')
        xlabel('time, (day)')
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        figure(2)
        subplot(2,3,1)
        plot(d4,grd.zt,'rs');
        ylim([0 6000])
        set(gca,'ydir','reverse','XAxisLocation','top')
        ylabel('depth (m)')
        xlabel('[^2^3^4Th_d] (mmol/L)')
        
        subplot(2,3,4)
        plot(d0,grd.zt,'rs');
        ylim([0 6000])
        set(gca,'ydir','reverse','XAxisLocation','top')
        ylabel('depth (m)')
        xlabel('[^2^3^0Th_d] (mmol/L)')
        
        subplot(2,3,2)
        plot(Th4,grd.zt,'rs');
        ylim([0 6000])
        set(gca,'ydir','reverse','XAxisLocation','top')
        %ylabel('depth (m)')
        xlabel('[^2^3^4Th_s] (\mumol/L)')
        
        subplot(2,3,5)
        plot(Th0,grd.zt,'rs');
        ylim([0 6000])
        set(gca,'ydir','reverse','XAxisLocation','top')
        %ylabel('depth (m)')
        xlabel('[^2^3^0Th_s] (\mumol/L)')
        
        subplot(2,3,3)
        plot(TH4,grd.zt,'rs');
        ylim([0 6000])
        set(gca,'ydir','reverse','XAxisLocation','top')
        %ylabel('depth (m)')
        xlabel('[^2^3^4Th_f] (\mumol/L)')
        
        subplot(2,3,6)
        plot(TH0,grd.zt,'rs');
        ylim([0 6000])
        set(gca,'ydir','reverse','XAxisLocation','top')
        xlabel('[^2^3^0Th_f] (\mumol/L)')
        
    end
    
end

