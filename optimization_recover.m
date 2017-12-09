load idealized_data.mat
load M3d.mat
load grid.mat
dpa = 365;
tmp = M3d;
tmp(:,:,4:end) = 0;
itop = find(tmp(:));

p.bs = 0.90;
p.bf = 0.75;

Ps0 = M3d*0;
Pf0 = M3d*0;
Chls0 = M3d*0;
Chlf0 = M3d*0;
Phs0 = M3d*0;
Phf0 = M3d*0;

Ps0(itop) = 16.30/60*4; 
Pf0(itop) = 27.98/60*4;
Chls0(itop) = 5.29/60*4;
Chlf0(itop) = 8.92/60*4;
Phs0(itop)  = 45.03/60*4;
Phf0(itop)  = 75.62/60*4;

p.Ps0 = Ps0(:);
p.Pf0 = Pf0(:);
p.Phs0 = Phs0(:);
p.Phf0 = Phf0(:);
p.Chls0 = Chls0(:);
p.Chlf0 = Chlf0(:);

p.Ps = idealized_data.Ps;
p.Pf = idealized_data.Pf;
p.Phs = idealized_data.Phs;
p.Phf = idealized_data.Phf;
p.Chls = idealized_data.Chls;
p.Chlf = idealized_data.Chlf;

p.Ps = p.Ps+0.5*randn(24,1);
p.Pf = p.Pf+0.5*randn(24,1);
p.Phs = p.Phs+0.5*randn(24,1);
p.Phf = p.Phf+0.5*randn(24,1);
p.Chls = p.Chls+0.5*randn(24,1);
p.Chlf = p.Chlf+0.5*randn(24,1);

p.Psstd = std(p.Ps);
p.Pfstd = std(p.Pf);
p.Phsstd = std(p.Phs);
p.Phfstd = std(p.Phf);
p.Chlsstd = std(p.Chls);
p.Chlfstd = std(p.Chlf);

dVt = grd.DZT3d.*grd.DXT3d.*grd.DYT3d;

p.k1 = 3/dpa;	
p.k2 = 150/dpa;	
p.d1 = 1/dpa;	
p.d2 = 1/dpa;

x0 = [p.k1;p.k2;p.d1;p.d2];
x0 = log(x0);
dx = sqrt(-1)*eps.^2*eye(length(x0));
%x0 = x0+dx(:,1);
L = @(x) neglogpost_recover(x,p,grd,M3d,dVt);

options = optimoptions(@fminunc,'Algorithm','trust-region',...
                       'GradObj','on','Hessian','on','Display','iter',...
                       'MaxFunEvals',5000,'MaxIter',5000,'TolX',1e-8,...
                       'DerivativeCheck','off','FinDiffType', ...
                       'central','TolFun',1e-8,'PrecondBandWidth',Inf);

[xhat,fval,exitflag] = fminunc(L,x0,options)

[f,dfdx,d2fdx2] = neglogpost_recover(xhat,p,grd,M3d,dVt)

sigmasquare = (2*f)/(2*138+1);
HH = d2fdx2./sigmasquare;
error = sqrt(diag(inv(HH)));
R.upbar = exp(xhat+error)-exp(xhat);
R.lowbar = exp(xhat)-exp(xhat-error);
R.xhat = exp(xhat)

fname = sprintf('Recover_results');
save(fname,'R');
