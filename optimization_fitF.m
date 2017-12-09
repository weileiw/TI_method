clc
clear all
Dab = 211; Dac = 1605; Dbc = 1394;
dpa = 365;

% =====================Reading data in.
filename='Finite_difference';
datasheet='49';

[ddata, dtext] = xlsread(filename,datasheet);
dname = dtext((2:length(dtext)),1);
dvalue = ddata(:,1);

for i = 1:length(dname)
    
    if strcmp(dname(i),'Chla')
        Chla = dvalue(i);
        
    elseif strcmp(dname(i),'ChlA')
        ChlA = dvalue(i);
        
    elseif strcmp(dname(i),'Pra')
        Pra = dvalue(i);
        
    elseif strcmp(dname(i),'PrA')
        PrA = dvalue(i);
        
    elseif strcmp(dname(i),'Chlb')
        Chlb = dvalue(i);
        
    elseif strcmp(dname(i),'ChlB')
        ChlB = dvalue(i);
        
    elseif strcmp(dname(i),'Prb')
        Prb = dvalue(i);
        
    elseif strcmp(dname(i),'PrB')
        PrB = dvalue(i);
        
    elseif strcmp(dname(i),'Chlc')
        Chlc = dvalue(i);
        
    elseif strcmp(dname(i),'ChlC')
        ChlC = dvalue(i);
        
    elseif strcmp(dname(i),'Prc')
        Prc = dvalue(i);
        
    elseif strcmp(dname(i),'PrC')
        PrC = dvalue(i);
        
    elseif strcmp(dname(i),'Pa')
        Pa = dvalue(i);
        
    elseif strcmp(dname(i),'PA')
        PA = dvalue(i);
        
    elseif strcmp(dname(i),'Pb')
        Pb = dvalue(i);
        
    elseif strcmp(dname(i),'PB')
        PB = dvalue(i);
        
    elseif strcmp(dname(i),'Pc')
        Pc = dvalue(i);
        
    elseif strcmp(dname(i),'PC')
        PC = dvalue(i);
        
    elseif strcmp(dname(i),'Fa')
        Fa = dvalue(i);
        
    elseif strcmp(dname(i),'FA')      
        FA = dvalue(i);
        
    elseif strcmp(dname(i),'Fb')
        Fb = dvalue(i);
        
    elseif strcmp(dname(i),'FB')
        FB = dvalue(i);
        
    elseif strcmp(dname(i),'Fc')
        Fc = dvalue(i);
        
    elseif strcmp(dname(i),'FC')
        FC = dvalue(i);
        
    elseif strcmp(dname(i), 'Fchla')
        Fchla=dvalue(i);
        
    elseif strcmp(dname(i), 'FchlA')
        FchlA=dvalue(i);
        
    elseif strcmp(dname(i), 'Fpra')
        Fpra=dvalue(i);
        
    elseif strcmp(dname(i), 'FprA')
        FprA=dvalue(i);
        
    elseif strcmp(dname(i), 'Fchlb')
        Fchlb=dvalue(i);
        
    elseif strcmp(dname(i), 'FchlB')
        FchlB=dvalue(i);
        
    elseif strcmp(dname(i), 'Fprb')
        Fprb=dvalue(i);
        
    elseif strcmp(dname(i), 'FprB')
        FprB=dvalue(i);
        
    elseif strcmp(dname(i), 'Fchlc')
        Fchlc=dvalue(i);
        
    elseif strcmp(dname(i), 'FchlC')
        FchlC=dvalue(i);
        
    elseif strcmp(dname(i), 'Fprc')
        Fprc=dvalue(i);
        
    elseif strcmp(dname(i), 'FprC')
        FprC=dvalue(i);
        
    end
    
end

za = 313; zb = 524; zc = 1918;

p.ChlA = 0.5*(zb*(ChlA+ChlB)+za*(ChlA-3*ChlB))/(zb-za);
p.ChlB = 0.5*(zc*(ChlB+ChlC)+zb*(ChlB-3*ChlC))/(zc-zb);
p.Chla = 0.5*(zb*(Chla+Chlb)+za*(Chla-3*Chlb))/(zb-za);
p.Chlb = 0.5*(zc*(Chlb+Chlc)+zb*(Chlb-3*Chlc))/(zc-zb);

p.PR1 = 0.5*(zb*(PrA+PrB)+za*(PrA-3*PrB))/(zb-za);
p.PR2 = 0.5*(zc*(PrB+PrC)+zb*(PrB-3*PrC))/(zc-zb);
p.pr1 = 0.5*(zb*(Pra+Prb)+za*(Pra-3*Prb))/(zb-za);
p.pr2 = 0.5*(zc*(Prb+Prc)+zb*(Prb-3*Prc))/(zc-zb);

p.P1 = 0.5*(zb*(PA+PB)+za*(PA-3*PB))/(zb-za);
p.P2 = 0.5*(zc*(PB+PC)+zb*(PB-3*PC))/(zc-zb);
p.p1 = 0.5*(zb*(Pa+Pb)+za*(Pa-3*Pb))/(zb-za);
p.p2 = 0.5*(zc*(Pb+Pc)+zb*(Pb-3*Pc))/(zc-zb);

p.Fla = (Fchlb-Fchla)/Dab;
p.Flb = (Fchlc-Fchlb)/Dbc;

p.Fra = (Fprb-Fpra)/Dab;
p.Frb = (Fprc-Fprb)/Dbc;

p.Fpa = (Fb-Fa)/Dab;
p.Fpb = (Fc-Fb)/Dbc;

p.FlA = (FchlB-FchlA)/Dab;
p.FlB = (FchlC-FchlB)/Dbc;

p.FrA = (FprB-FprA)/Dab;
p.FrB = (FprC-FprB)/Dbc;

p.FpA = (FB-FA)/Dab;
p.FpB = (FC-FB)/Dbc;

p.k1 = 3/dpa; p.k2 = 150/dpa;
p.d1 = 1/dpa; p.d2 = 1/dpa;
x0 = [p.k1;p.k2;p.d1;p.d2];
x0 = log(x0);
%load xhat_sedi.mat
x0 = [-23.7339;-3.4405;-5.5527;-5.5351];
L = @(x) neglogpost_fitF(x,p);

options = optimoptions(@fminunc,'Algorithm','trust-region',...
                       'GradObj','on','Hessian','on','Display','iter',...
                       'MaxFunEvals',5000,'MaxIter',5000,'TolX',1e-8,...
                      'DerivativeCheck','off','FinDiffType', ...
                      'central','TolFun',1e-8,'PrecondBandWidth',Inf);

[xhat,fval,exitflag] = fminunc(L,x0,options)

[f,dfdx,d2fdx2] = neglogpost_fitF(xhat,p);

s2 = 2*f/(2*12+1);
HH = d2fdx2./s2;
error = sqrt(diag(inv(HH)));
R.upbar = exp(xhat+error)-exp(xhat);
R.lowbar = exp(xhat)-exp(xhat-error);
R.xhat = exp(xhat);

% $$$ fname = sprintf('xhat_sedi_v3');
% $$$ save(fname,'R');

