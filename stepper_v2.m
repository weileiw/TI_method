clear all
close all
load M3d.mat
load grid.mat
p.bs = 0.90;
p.bf = 0.75;
figure(1)
dpa = 365; 

k1 = 3/dpa;	
k2 = 150/dpa;	
d1 = 1/dpa;	
d2 = 1/dpa;	
d3 = 1/dpa
p.k1 = k1; p.k2 = k2; 
p.d1 = d1; p.d2 = d2;
p.d3 = d3;
Ps   = 0*ones(24,1);
Pf   = 0*ones(24,1);
Chls = 0*ones(24,1); 	
Chlf = 0*ones(24,1);	
Phs  = 0*ones(24,1);	
Phf  = 0*ones(24,1);

tmp = M3d;
tmp(:,:,4:end) = 0;
itop = find(tmp(:));

Ps0 = M3d*0;
Pf0 = M3d*0;
Chls0 = M3d*0;
Chlf0 = M3d*0;
Phs0 = M3d*0;
Phf0 = M3d*0;
% flux estimated from sediment traps divived
% by eqphotic zone depth.
Ps0(itop) = 16.30/60*4; 
Pf0(itop) = 27.98/60*4;
Chls0(itop) = 5.29/60*4;
Chlf0(itop) = 8.92/60*4;
Phs0(itop)  = 45.03/60*4;
Phf0(itop)  = 75.62/60*4;

PFdiv_s = buildPFD_s(M3d,p,grd);
PFdiv_c = buildPFD_Chl(M3d,p,grd);
PFdiv_f = buildPFD_f(M3d,p,grd);
PFdiv_s(24,24) = PFdiv_s(23,23);
PFdiv_c(24,24) = PFdiv_c(23,23);
PFdiv_f(24,24) = PFdiv_f(23,23);

t = 0;	
dt = 1/24;
I = speye(length(Ps));
A = I+(dt/2)*PFdiv_s;
B = I-(dt/2)*PFdiv_s;
C = I+(dt/2)*PFdiv_f;
D = I-(dt/2)*PFdiv_f;
E = I+(dt/2)*PFdiv_c;
F = I-(dt/2)*PFdiv_c;

FA = mfactor(A);
FC = mfactor(C);
FE = mfactor(E);

for ik = 1:24*2*365
    
    dPsdt   = Pf*k2-Ps*k1-Ps*d2+Ps0(:);
    dChlsdt = Chlf*k2-Chls*k1-Chls*d1+Chls0(:);
    dPhsdt  = Phf*k2-Phs*k1-Phs*d3+Chls*d1+Phs0(:);
    
    dPfdt   = Ps*k1-Pf*k2+Pf0(:);
    dChlfdt = Chls*k1-Chlf*k2+Chlf0(:);
    dPhfdt  = Phs*k1-Phf*k2+Phf0(:);
    
    Ps   = mfactor(FA,(B*Ps+dPsdt*dt));
    Chls = mfactor(FE,(F*Chls+dChlsdt*dt));
    Phs  = mfactor(FA,(B*Phs+dPhsdt*dt));
    Pf   = mfactor(FC,(D*Pf+dPfdt*dt));
    Chlf = mfactor(FC,(D*Chlf+dChlfdt*dt));
    Phf  = mfactor(FC,(D*Phf+dPhfdt*dt));
    t = t+dt;
    
    if mod(ik,48) == 0
        
        set(0,'CurrentFigure',1);
        subplot(3,2,1)
        plot(t,Ps(24),'rs'); hold on; drawnow	
        ylabel('[Particle-slow] (mmol/L)')
        xlabel('time, (day)')
        
        subplot(3,2,2)
        plot(t,Pf(24),'rs');hold on; drawnow
        ylabel('[Particle-fast] (mmol/L)')
        xlabel('time, (day)')
        
        subplot(3,2,3)
        plot(t,Chls(24),'rs');hold on; drawnow
        ylabel('[Chla-slow] (\mu mol/L)')
        xlabel('time, (day)')
        
        subplot(3,2,4)
        plot(t,Chlf(24),'rs'); hold on;drawnow
        ylabel('[Chla-fast] (\mu mol/L)')
        xlabel('time, (day)')
        
        subplot(3,2,5)
        plot(t,Phs(24),'rs');hold on;
        ylabel('[Pheo-slow] (\mu mol/L)')
        xlabel('time, (day)')
        
        subplot(3,2,6)
        plot(t,Phf(24),'rs');hold on;drawnow
        ylabel('[Pheo-fast] (\mu mol/L)')
        xlabel('time, (day)')
        
    end
    
end

figure(2)
%plot(Chls,grd.dzt)
subplot(2,3,1)
plot(Ps,grd.zt,'rs');
ylim([0 6000])
set(gca,'ydir','reverse','XAxisLocation','top')
ylabel('depth (m)')
xlabel('[Particle-slow] (mmol/L)')

subplot(2,3,4)
plot(Pf,grd.zt,'rs');
ylim([0 6000])
set(gca,'ydir','reverse','XAxisLocation','top')
ylabel('depth (m)')
xlabel('[Particle-fast] (mmol/L)')

subplot(2,3,2)
plot(Chls,grd.zt,'rs');
ylim([0 6000])
set(gca,'ydir','reverse','XAxisLocation','top')
%ylabel('depth (m)')
xlabel('[Chla-slow] (\mu mol/L)')

subplot(2,3,5)
plot(Chlf,grd.zt,'rs');
ylim([0 6000])
set(gca,'ydir','reverse','XAxisLocation','top')
%ylabel('depth (m)')
xlabel('[Chla-fast] (\mu mol/L)')

subplot(2,3,3)
plot(Phs,grd.zt,'rs');
ylim([0 6000])
set(gca,'ydir','reverse','XAxisLocation','top')
%ylabel('depth (m)')
xlabel('[Pheo-slow] (\mu mol/L)')

subplot(2,3,6)
plot(Phf,grd.zt,'rs');
ylim([0 6000])
set(gca,'ydir','reverse','XAxisLocation','top')
%ylabel('depth (m)')
xlabel('[Pheo-fast] (\mu mol/L)')
 
idealized_data.Ps = Ps;
idealized_data.Pf = Pf;
idealized_data.Phs = Phs;
idealized_data.Phf = Phf;
idealized_data.Chls = Chls;
idealized_data.Chlf = Chlf;

fname = sprintf('idealized_data');
save(fname,'idealized_data')



