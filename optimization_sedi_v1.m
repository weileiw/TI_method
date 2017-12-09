Dab = 211; Dac = 1605; Dbc = 1394;
dpa = 365;
spd = 24*60*60;
spa = dpa*spd;
% ===========================Reading parameters in, the limit of each
% paremeters were initially set in a arbitray way, and were changed
% accordingly based on the first run results.


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

ChlA = 0.5*(zb*(ChlA+ChlB)+za*(ChlB-3*ChlA))/(zb-za);
ChlB = 0.5*(zc*(ChlB+ChlC)+zb*(ChlC-3*ChlB))/(zc-zb);
Chla = 0.5*(zb*(Chla+Chlb)+za*(Chlb-3*Chla))/(zb-za);
Chlb = 0.5*(zc*(Chlb+Chlc)+zb*(Chlc-3*Chlb))/(zc-zb);

PR1 = 0.5*(zb*(PrA+PrB)+za*(PrB-3*PrA))/(zb-za);
PR2 = 0.5*(zc*(PrB+PrC)+zb*(PrC-3*PrB))/(zc-zb);
pr1 = 0.5*(zb*(Pra+Prb)+za*(Prb-3*Pra))/(zb-za);
pr2 = 0.5*(zc*(Prb+Prc)+zb*(Prc-3*Prb))/(zc-zb);

P1 = 0.5*(zb*(PA+PB)+za*(PB-3*PA))/(zb-za);
P2 = 0.5*(zc*(PB+PC)+zb*(PC-3*PB))/(zc-zb);
p1 = 0.5*(zb*(Pa+Pb)+za*(Pb-3*Pa))/(zb-za);
p2 = 0.5*(zc*(Pb+Pc)+zb*(Pc-3*Pb))/(zc-zb);

Ps = @(z) log(z/5899.2)./(-0.612);
Pf = @(z) log(grd.zt./5466)./(-21.7);
Phs = @(z) log(grd.zt./3953.9)./(-0.19);
Phf = @(z) log(grd.zt./3066.7)./(-5.313);
Chls = @(z) log(grd.zt./4270)./(-1.616);
Chlf = @(z) log(grd.zt./5878.3)./(-60.52);

p.ChlA = ChlA/std([ChlA,ChlB]);
p.ChlB = ChlB/std([ChlA,ChlB]);
p.Chla = Chla/std([Chla,Chlb]);
p.Chlb = Chlb/std([Chla,Chlb]);

p.PR1 = PR1/std([PR1,PR2]);
p.PR2 = PR2/std([PR1,PR2]);
p.pr1 = pr1/std([pr1,pr2]);
p.pr2 = pr2/std([pr1,pr2]);

P1 = P1/std([P1,P2]);
P2 = P2/std([P1,P2]);
p1 = p1/std([p1,p2]);
p2 = p1/std([p1,p2]);

Fla = (Fchla-Fchlb)/Dab;
Flb = (Fchlb-Fchlc)/Dbc;

Fra = (Fpra-Fprb)/Dab;
Frb = (Fprb-Fprc)/Dbc;

Fpa = (Fa-Fb)/Dab;
Fpb = (Fb-Fc)/Dbc;

FlA = (FchlA-FchlB)/Dab;
FlB = (FchlB-FchlC)/Dbc;

FrA = (FprA-FprB)/Dab;
FrB = (FprB-FprC)/Dbc;

FpA = (FA-FB)/Dab;
FpB = (FB-FC)/Dbc;

p.Fla = Fla/std([Fla,Flb]);
p.Flb = Flb/std([Fla,Flb]);

p.Fra = Fra/std([Fra,Frb]);
p.Frb = Frb/std([Fra,Frb]);

p.Fpa = Fpa/std([Fpa,Fpb]);
p.Fpb = Fpb/std([Fpa,Fpb]);

p.FlA = FlA/std([FlA,FlB]);
p.FlB = FlB/std([FlA,FlB]);

p.FrA = FrA/std([FrA,FrB]);
p.FrB = FrB/std([FrA,FrB]);

p.FpA = FpA/std([FpA,FpB]);
p.FpB = FpB/std([FpA,FpB]);

x0 = [3;150;10;10];
x0 = [log(x0)];
%dx = sqrt(-1)*speye(4)*eps.^2;
%x0 = x0+dx(:,4);
keyboard
L = @(x) neglogpost_sedi_v1(x,p);

options = optimoptions(@fminunc,'Algorithm','trust-region',...
                       'GradObj','on','Hessian','on','Display','iter',...
                       'MaxFunEvals',5000,'MaxIter',5000,'TolX',1e-7,...
                      'DerivativeCheck','off','FinDiffType', ...
                      'central','TolFun',1e-7,'PrecondBandWidth',Inf);

[xhat,fval,exitflag] = fminunc(L,x0,options)
%options = optimset('Display','iter','TolFun',1e-6,'MaxIter',5000,'MaxFunEvals',5000);
%[xhat,fval,exitflag] = fminsearch(L,x0,options)

[f,dfdx,d2fdx2] = neglogpost_sedi_v1(xhat,p);

s2 = (2*f)/(2*12+1);
HH = d2fdx2./s2;
error = sqrt(diag(inv(HH)));
R.upbar = exp(xhat+error)-exp(xhat);
R.lowbar = exp(xhat)-exp(xhat-error);
R.xhat = exp(xhat);

fname = sprintf('xhat_sedi_v1');
save(fname,'R');
