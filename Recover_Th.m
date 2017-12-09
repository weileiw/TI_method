load Steady_Th234_Th230.mat
load M3d.mat
load grid.mat
dpa = 365;
tmp = M3d;

p.Th0   = Th0;
p.Th4   = Th4;
p.TH0   = TH0;
p.TH4   = TH4
p.d230  = d230;
p.d234  = d234;

p.Th0_n   = Th0.*(1+0.25*randn(24,1));
p.Th4_n   = Th4.*(1+0.25*randn(24,1));
p.TH0_n   = TH0.*(1+0.25*randn(24,1));
p.TH4_n   = TH4.*(1+0.25*randn(24,1));
p.d230_n  = d230.*(1+0.25*randn(24,1));
p.d234_n  = d234.*(1+0.25*randn(24,1));

p.Th0_std  = std(p.Th0_n);
p.Th4_std  = std(p.Th4_n);
p.TH0_std  = std(p.TH0_n);
p.TH4_std  = std(p.TH4_n);
p.d230_std = std(p.d230_n);
p.d234_std = std(p.d234_n);

dVt = grd.DZT3d.*grd.DXT3d.*grd.DYT3d;

p.aggregation      = 3/dpa;	
p.disagregation    = 150/dpa;	
p.remineralization = 1/dpa;	
p.adsorption       = 1/dpa;
p.desorption       = 2/dpa;

x0 = [p.aggregation;p.disagregation;p.remineralization;...
      p.adsorption;p.desorption];
x0 = log(x0);
dx = sqrt(-1)*eps.^2*eye(length(x0));
%x0 = x0+dx(:,1);
L = @(x) EQ_Th(x,p,grd,M3d,dVt);

options = optimoptions(@fminunc,'Algorithm','trust-region',...
                       'GradObj','on','Hessian','on','Display','iter',...
                       'MaxFunEvals',5000,'MaxIter',5000,'TolX',1e-8,...
                       'DerivativeCheck','off','FinDiffType', ...
                       'central','TolFun',1e-8,'PrecondBandWidth',Inf);

[xhat,fval,exitflag] = fminunc(L,x0,options)

[f,dfdx,d2fdx2] = EQ_Th(xhat,p,grd,M3d,dVt)

sigmasquare = (2*f)/(2*138+1);
HH = d2fdx2./sigmasquare;
error = sqrt(diag(inv(HH)));
R.upbar = exp(xhat+error)-exp(xhat);
R.lowbar = exp(xhat)-exp(xhat-error);
R.xhat = exp(xhat)

fname = sprintf('Recover_results');
save(fname,'R');
