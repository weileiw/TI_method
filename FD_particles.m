% this script intends to generate partile distribution by using 
% finite difference method.
clear all
close all
load M3d.mat
load grid.mat
p.bf = 0.95;

dpa        = 365;     % day per year

p.w        = 150*dpa;     % m y^-1;
beta_nega1 = 1/dpa;       % Remi y^-1;
beta2      = 3/dpa;       % Aggre y^-1;
beta_nega2 = 150/dpa;     % disaggre y^-1;
p.k2       = beta_nega2;

% depth 
z = linspace(110,4900,100);
% large particle concentration at 110m 
PFD = buildPFD_cons_SV(M3d,p,grd);
PFD(end,end) = PFD(end-1,end-1);
%P_l = zeros(length(z),1);
P_l = zeros(24,1);
P_s = P_l;
P_l(1) = 1e-6;

I  = speye(length(P_l));

t  = 0;
dt = 1/24;
nstep = 10*365/dt;

A = I+(dt/2)*PFD;
B = I-(dt/2)*PFD;
FA = mfactor(A);

figure(1)
for jl = 1:nstep

    dP_sdt = beta_nega2*P_l - (beta_nega1+beta2)*P_s;
    dP_ldt = beta2*P_s - beta_nega2*P_l;

    P_l = mfactor(FA,(B*[1e-6;P_l(2:end)]+dP_ldt*dt));
    P_s = P_s+dP_sdt*dt;
    t = t+dt;
    if mod(jl,30) == 0
        set(0,'CurrentFigure',1);
        subplot(2,1,1)
        plot(t,P_l(20),'rs'); hold on; drawnow  

        subplot(2,1,2)
        plot(t,P_s(20),'rs'); hold on; drawnow  
        
        figure(2)
        subplot(1,2,1)
        plot(P_l,grd.zt,'rs');
        ylim([0 5000])
        set(gca,'ydir','reverse','XAxisLocation','top')
        ylabel('depth (m)')
        xlabel('large particle concentration')
        subplot(1,2,2)
        plot(P_s,grd.zt,'rs');
        ylim([0 5000])
        set(gca,'ydir','reverse','XAxisLocation','top')
        ylabel('depth (m)')
        xlabel('small particle concentration')

    end
end
