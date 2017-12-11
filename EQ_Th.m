function [f,dfdx,d2fdx2] = EQ_Th(x,p,grd,M3d,dVt)

nip = length(x);
%dx = eps^2*sqrt(-1)*speye(nip);
%x = x+dx(:,5);
p.aggregation      = exp(x(1));
p.disagregation    = exp(x(2));
p.remineralization = exp(x(3));
p.adsorption       = exp(x(4));
p.desorption       = exp(x(5));

[M,D,D2] = Th_cycle(p,grd,M3d);
idata = 24;
d234 = M(1:idata);
Th4  = M(1*idata+1:2*idata);
TH4  = M(2*idata+1:3*idata);
d230 = M(3*idata+1:4*idata);
Th0  = M(4*idata+1:5*idata);
TH0  = M(5*idata+1:6*idata);
P_s  = M(6*idata+1:7*idata);
P_l  = M(7*idata+1:8*idata);

dVt = dVt(1:24);

W = repmat(dVt(:),8,1);
W = d0(W);

%W = d0(dVt(:)./sum(dVt(:)));

e1 = (d234-p.d234_n)/p.d234_std;
e2 = (Th4 -p.Th4_n)/p.Th4_std;
e3 = (TH4 -p.TH4_n)/p.TH4_std;
e4 = (d230-p.d230_n)/p.d230_std;
e5 = (Th0 -p.Th0_n)/p.Th0_std;
e6 = (TH0 -p.TH0_n)/p.TH0_std;
e7 = (P_s -p.P_s_n)/p.P_s_std;
e8 = (P_l -p.P_l_n)/p.P_l_std;


e = [e1;e2;e3;e4;e5;e6;e7;e8];

%f = 0.5*(e1.'*W*e1+e2.'*W*e2+e3.'*W*e3+e4.'*W*e4+e5.'*W*e5+e6.'*W* ...
%         e6);%+0.00011*(ep.'*Wp*ep);
f = 0.5*e.'*W*e;
D = [D(1:idata,:)/p.d234_std;...
     D(1*idata+1:2*idata,:)/p.Th4_std;...
     D(2*idata+1:3*idata,:)/p.TH4_std;...
     D(3*idata+1:4*idata,:)/p.d230_std;...
     D(4*idata+1:5*idata,:)/p.Th0_std;...
     D(5*idata+1:6*idata,:)/p.TH0_std;...
     D(6*idata+1:7*idata,:)/p.P_s_std;...
     D(7*idata+1:8*idata,:)/p.P_l_std];

dpdx = diag(exp(x));
dedx = D*dpdx;
dfdx = e.'*W*dedx;

D2 = [D2(1:idata,:)/p.d234_std;...
      D2(1*idata+1:2*idata,:)/p.Th4_std;...
      D2(2*idata+1:3*idata,:)/p.TH4_std;...
      D2(3*idata+1:4*idata,:)/p.d230_std;...
      D2(4*idata+1:5*idata,:)/p.Th0_std;...
      D2(5*idata+1:6*idata,:)/p.TH0_std;...
      D2(6*idata+1:7*idata,:)/p.P_s_std;...
      D2(7*idata+1:8*idata,:)/p.P_l_std];

H(1,:) = e.'*W*(D2(:,1:5));
H(2,:) = e.'*W*(D2(:,6:10));
H(3,:) = e.'*W*(D2(:,11:15));
H(4,:) = e.'*W*(D2(:,16:20));
H(5,:) = e.'*W*(D2(:,21:25));

H = 0.5*dpdx*(H+H.')*dpdx+diag(e.'*W*dedx);
d2fdx2 = dedx.'*W*dedx+H;

function [M,D,D2] = Th_cycle(p,grd,M3d)

k1 = p.aggregation;
k2 = p.disagregation;
r  = p.remineralization;
a  = p.adsorption;
d  = p.desorption;
n234  = p.n234;
n230  = p.n230;
U238  = p.U238;
U234  = p.U234;

PFD = buildPFD_cons_SV(M3d,p,grd);
PFD(24,24) = PFD(23,23);

I = speye(length(U234));
INIT_4 = [ 2.5;0*U238(2:end)];
INIT_0 = [1e-3;0*U234(2:end)];
INIT_P = [1e-6;0*U234(2:end)];
rhs = [U238*n234;0*U238;PFD*INIT_4;U234*n230;0*U234;PFD*INIT_0;0*U234;PFD*INIT_P];

J = [[ (n234+a)*I,       -(d+r)*I,             0*I, 0*I,0*I,0*I,0*I,0*I];...
     [       -a*I,(d+r+k1+n234)*I,           -k2*I, 0*I,0*I,0*I,0*I,0*I];...
     [        0*I,          -k1*I, (k2+n234)*I+PFD, 0*I,0*I,0*I,0*I,0*I];...
     [0*I,0*I,0*I, (n230+a)*I,        -(d+r)*I,             0*I,0*I,0*I];...
     [0*I,0*I,0*I,       -a*I, (d+r+k1+n230)*I,           -k2*I,0*I,0*I];...
     [0*I,0*I,0*I,        0*I,           -k1*I, (k2+n230)*I+PFD,0*I,0*I];...
     [0*I,0*I,0*I,0*I,0*I,0*I, (r+k1)*I,     -k2*I];...
     [0*I,0*I,0*I,0*I,0*I,0*I,    -k1*I,  k2*I+PFD]];

FJ = mfactor(J);
M = mfactor(FJ,rhs);

dJdk1 = [[0*I, 0*I, 0*I, 0*I, 0*I, 0*I, 0*I, 0*I];...
         [0*I,   I, 0*I, 0*I, 0*I, 0*I, 0*I, 0*I];...
         [0*I,  -I, 0*I, 0*I, 0*I, 0*I, 0*I, 0*I];...	
         [0*I, 0*I, 0*I, 0*I, 0*I, 0*I, 0*I, 0*I];...
         [0*I, 0*I, 0*I, 0*I,   I, 0*I, 0*I, 0*I];...
         [0*I, 0*I, 0*I, 0*I,  -I, 0*I, 0*I, 0*I];...
         [0*I, 0*I, 0*I, 0*I, 0*I, 0*I,   I, 0*I];...
         [0*I, 0*I, 0*I, 0*I, 0*I, 0*I,  -I, 0*I]];

dJdk2 = [[0*I, 0*I, 0*I, 0*I, 0*I, 0*I, 0*I, 0*I];...
         [0*I, 0*I,  -I, 0*I, 0*I, 0*I, 0*I, 0*I];...
         [0*I, 0*I,   I, 0*I, 0*I, 0*I, 0*I, 0*I];...	
         [0*I, 0*I, 0*I, 0*I, 0*I, 0*I, 0*I, 0*I];...
         [0*I, 0*I, 0*I, 0*I, 0*I,  -I, 0*I, 0*I];...
         [0*I, 0*I, 0*I, 0*I, 0*I,   I, 0*I, 0*I];...
         [0*I, 0*I, 0*I, 0*I, 0*I, 0*I, 0*I,  -I];...
         [0*I, 0*I, 0*I, 0*I, 0*I, 0*I, 0*I,   I]];

dJda = [[  I, 0*I, 0*I, 0*I, 0*I, 0*I, 0*I, 0*I];...
        [ -I, 0*I, 0*I, 0*I, 0*I, 0*I, 0*I, 0*I];...
        [0*I, 0*I, 0*I, 0*I, 0*I, 0*I, 0*I, 0*I];...	
        [0*I, 0*I, 0*I,   I, 0*I, 0*I, 0*I, 0*I];...
        [0*I, 0*I, 0*I,  -I, 0*I, 0*I, 0*I, 0*I];...
        [0*I, 0*I, 0*I, 0*I, 0*I, 0*I, 0*I, 0*I];...
        [0*I, 0*I, 0*I, 0*I, 0*I, 0*I, 0*I, 0*I];...
        [0*I, 0*I, 0*I, 0*I, 0*I, 0*I, 0*I, 0*I]];

dJdd = [[0*I,  -I, 0*I, 0*I, 0*I, 0*I, 0*I, 0*I];...
        [0*I,   I, 0*I, 0*I, 0*I, 0*I, 0*I, 0*I];...
        [0*I, 0*I, 0*I, 0*I, 0*I, 0*I, 0*I, 0*I];...	
        [0*I, 0*I, 0*I, 0*I,  -I, 0*I, 0*I, 0*I];...
        [0*I, 0*I, 0*I, 0*I,   I, 0*I, 0*I, 0*I];...
        [0*I, 0*I, 0*I, 0*I, 0*I, 0*I, 0*I, 0*I];...
        [0*I, 0*I, 0*I, 0*I, 0*I, 0*I, 0*I, 0*I];...
        [0*I, 0*I, 0*I, 0*I, 0*I, 0*I, 0*I, 0*I]];

dJdr = [[0*I,  -I, 0*I, 0*I, 0*I, 0*I, 0*I, 0*I];...
        [0*I,   I, 0*I, 0*I, 0*I, 0*I, 0*I, 0*I];...
        [0*I, 0*I, 0*I, 0*I, 0*I, 0*I, 0*I, 0*I];...	
        [0*I, 0*I, 0*I, 0*I,  -I, 0*I, 0*I, 0*I];...
        [0*I, 0*I, 0*I, 0*I,   I, 0*I, 0*I, 0*I];...
        [0*I, 0*I, 0*I, 0*I, 0*I, 0*I, 0*I, 0*I];...
        [0*I, 0*I, 0*I, 0*I, 0*I, 0*I,   I, 0*I];...
        [0*I, 0*I, 0*I, 0*I, 0*I, 0*I, 0*I, 0*I]];


D = mfactor(FJ,-[dJdk1*M,...
	         dJdk2*M,...
		 dJdr*M,...
		 dJda*M,...
                 dJdd*M]);

D2 = mfactor(FJ,[-2*dJdk1*D,...
	         -2*dJdk2*D,...
		 -2*dJdr*D,...
		 -2*dJda*D,...
                 -2*dJdd*D]);

