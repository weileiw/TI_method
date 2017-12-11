clear all
close all
addpath('~/MOCM/WEILEI/myfunc')
load M3d.mat
load grid.mat
load Steady_Th234_Th230.mat
p.bf = 0.95;
%figure(1)
dpa = 365; 
n234 = 10.5/dpa; 
n230 = 9.19e-6/dpa;

U238   = 2400*ones(24,1);
U234   = 2760*ones(24,1);
P_l    = zeros(24,1);
P_s    = P_l;
%d230  = 0*ones(24,1);
%d234  = 0*ones(24,1);
%TH4 = 0*ones(24,1); 	
%TH0 = 0*ones(24,1);
%Th4  = 0*ones(24,1);	
%Th0  = 0*ones(24,1);

k1 = 3/dpa;   % aggregation
k2 = 150/dpa; % disagregation	
a  = 0.5/dpa; % adsorption
d  = 2/dpa;   % desorption
r  = 1/dpa;   % remineralization
              %p.k2 = k2;
PFD = buildPFD_cons_SV(M3d,p,grd);
PFD(end,end) = PFD(end-1,end-1);

t  = 0;	
dt = 1/24;
I  = speye(length(U238));
A  = I+(dt/2)*PFD;
B  = I-(dt/2)*PFD;

FA = mfactor(A);

nstep = 50*365/dt;
for ik = 1:nstep
    
    dP_sdt  = P_l*k2 - P_s*(k1+r);
    dP_ldt  = P_s*k1 - P_l*k2;%+[1e-6;0*P_s(2:end)];  
    
    dd234dt = (U238-d234)*n234-d234*a+Th4*(d+r);
    dTh4dt  = d234*a+TH4*k2-Th4*(k1+r+d+n234);
    dTH4dt  = Th4*k1-TH4*(k2+n234);%+[2.5; 0*TH4(2:end)];
    
    dd230dt = (U234-d230)*n230-d230*a+Th0*(d+r); 
    dTh0dt  = d230*a+TH0*k2-Th0*(k1+r+d+n230);
    dTH0dt  = Th0*k1-TH0*(k2+n230);%+[0.001;0*TH0(2:end)];
    
    P_s  = P_s+dP_sdt*dt;
    d234 = d234+dd234dt*dt;
    d230 = d230+dd230dt*dt;
    Th4  = Th4+dTh4dt*dt;
    Th0  = Th0+dTh0dt*dt;

    P_l  = mfactor(FA,(B*[1e-6;P_l(2:end)]+dP_ldt*dt));
    TH4  = mfactor(FA,(B*[ 2.5;TH4(2:end)]+dTH4dt*dt));
    TH0  = mfactor(FA,(B*[1e-3;TH0(2:end)]+dTH0dt*dt));

    t = t+dt;
    if mod(ik,24*50) == 0
        fprintf('model time is %d days; \n', ik/24);
    end
    %if mod(ik,48) == 0
    
    %set(0,'CurrentFigure',1);
    %subplot(3,2,1)
    %plot(t,d234(24),'rs'); hold on; drawnow
    %ylabel('[^2^3^4Th_d] (dpm/m^3)')
    %xlabel('time, (day)')
    
    %subplot(3,2,2)
    %plot(t,d230(24),'rs');hold on; drawnow
    %ylabel('[^2^3^0Th_d (dpm/m^3)')
    %xlabel('time, (day)')
    
    %subplot(3,2,3)
    %plot(t,TH4(24),'rs');hold on; drawnow
    %ylabel('[^2^3^4Th_f] (dpm/m^3)')
    %xlabel('time, (day)')
    
    %subplot(3,2,4)
    %plot(t,TH0(24),'rs'); hold on;drawnow
    %ylabel('[^2^3^0Th_f] (dpm/m^3)')
    %xlabel('time, (day)')
    
    %subplot(3,2,5)
    %plot(t,Th4(24),'rs');hold on;
    %ylabel('[^2^3^4Th_s] (dpm/m^3)')
    %xlabel('time, (day)')
    
    %subplot(3,2,6)
    %plot(t,Th0(24),'rs');hold on;drawnow
    %ylabel('[^2^3^0Th_s] (dpm/m^3)')
    %xlabel('time, (day)')
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %figure(2)
    %subplot(2,3,1)
    %plot(d234,grd.zt,'rs');
    %ylim([0 6000])
    %set(gca,'ydir','reverse','XAxisLocation','top')
    %ylabel('depth (m)')
    %xlabel('[^2^3^4Th_d] (mmol/L)')
    
    %subplot(2,3,4)
    %plot(d230,grd.zt,'rs');
    %ylim([0 6000])
    %set(gca,'ydir','reverse','XAxisLocation','top')
    %ylabel('depth (m)')
    %xlabel('[^2^3^0Th_d] (mmol/L)')
    
    %subplot(2,3,2)
    %plot(Th4,grd.zt,'rs');
    %ylim([0 6000])
    %set(gca,'ydir','reverse','XAxisLocation','top')
				%%ylabel('depth (m)')
    %xlabel('[^2^3^4Th_s] (\mumol/L)')
    
    %subplot(2,3,5)
    %plot(Th0,grd.zt,'rs');
    %ylim([0 6000])
    %set(gca,'ydir','reverse','XAxisLocation','top')
				%%ylabel('depth (m)')
    %xlabel('[^2^3^0Th_s] (\mumol/L)')
    
    %subplot(2,3,3)
    %plot(TH4,grd.zt,'rs');
    %ylim([0 6000])
    %set(gca,'ydir','reverse','XAxisLocation','top')
				%%ylabel('depth (m)')
    %xlabel('[^2^3^4Th_f] (\mumol/L)')
    
    %subplot(2,3,6)
    %plot(TH0,grd.zt,'rs');
    %ylim([0 6000])
    %set(gca,'ydir','reverse','XAxisLocation','top')
    %xlabel('[^2^3^0Th_f] (\mumol/L)')
    
  %end
  
end

fname = sprintf('Steady_Th234_Th230');
save(fname,'Th0','Th4','TH0','TH4','d230','d234','P_l','P_s')