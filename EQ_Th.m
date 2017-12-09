function [f,dfdx,d2fdx2] = EQ_Th(x,p,grd,M3d,dVt)

nip = length(x);

p.aggregation      = exp(x(1));
p.disagregation    = exp(x(2));
p.remineralization = exp(x(3));
p.adsorption       = exp(x(4));
p.desorption       = exp(x(5));

[M,D,D2] = Th_cycle(p,grd,M3d);

d234 = M(1:24);
Th4  = M(25:48);
TH4  = M(49:72);
d230 = M(73:96);
Th0  = M(97:120);
TH0  = M(121:144);

dVt = dVt(1:24);

W = repmat(dVt(:),6,1);
W = d0(W);

%W = d0(dVt(:)./sum(dVt(:)));

e1 = (d234-p.d234)/p.d234_std;
e2 = (Th4 -p.Th4)/p.Th4_std;
e3 = (TH4 -p.TH4)/p.TH4_std;
e4 = (d230-p.d230)/p.d230_std;
e5 = (Th0 -p.Th0)/p.Th0_std;
e6 = (TH0 -p.TH0)/p.TH0_std;

e = [e1;e2;e3;e4;e5;e6];

%f = 0.5*(e1.'*W*e1+e2.'*W*e2+e3.'*W*e3+e4.'*W*e4+e5.'*W*e5+e6.'*W* ...
%         e6);%+0.00011*(ep.'*Wp*ep);
f = 0.5*e.'*W*e;
D = [D(1:24,:)/p.d234_std;...
     D(25:48,:)/p.Th4_std;...
     D(49:72,:)/p.TH4_std;...
     D(73:96,:)/p.d230_std;...
     D(97:120,:)/p.Th0_std;...
     D(121:144,:)/p.TH0_std];

dpdx = diag(exp(x));
dedx = D*dpdx;
dfdx = e.'*W*dedx;

D2 = [D2(1:24,:)/p.d234_std;...
      D2(25:48,:)/p.Th4_std;...
      D2(49:72,:)/p.TH4_std;...
      D2(73:96,:)/p.d230_std;...
      D2(97:120,:)/p.Th0_std;...
      D2(121:144,:)/p.TH0_std];

H(1,:) = e.'*W*(D2(:,1:5));
H(2,:) = e.'*W*(D2(:,6:10));
H(3,:) = e.'*W*(D2(:,11:15));
H(4,:) = e.'*W*(D2(:,16:20));
H(5,:) = e.'*W*(D2(:,21:25));

H = 0.5*dpdx*(H+H.')*dpdx+diag(e.'*W*dedx);
d2fdx2 = dedx.'*W*dedx+H;
keyboard
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

I = speye(length(d234));
INIT_4 = [  2.5,0*U238(2:end)];
INIT_0 = [0.001,0*U234(2:end)];

rhs = [U238*n234;0*U238;INIT_4;U234*n230;0*U234;INIT_0];

J = [[ (n234+a)*I,       -(d+r)*I,             0*I, 0*I,0*I,0*I];...
     [       -a*I,(d+r+k1+n234)*I,           -k2*I, 0*I,0*I,0*I];...
     [        0*I,          -k1*I, (k2+n234)*I+PFD, 0*I,0*I,0*I];...
     [0*I,0*I,0*I, (n230+a)*I,        -(d+r)*I,            0*I];...
     [0*I,0*I,0*I,       -a*I, (d+r+k1+n230)*I,          -k2*I];...
     [0*I,0*I,0*I,        0*I,           -k1*I,(k2+n230)*I+PFD]];

FJ = mfactor(J);
M = mfactor(FJ,rhs);

dJdk1 = [[0*I, 0*I, 0*I, 0*I, 0*I, 0*I];...
         [0*I,   I, 0*I, 0*I, 0*I, 0*I];...
         [0*I,  -I, 0*I, 0*I, 0*I, 0*I];...	
         [0*I, 0*I, 0*I, 0*I, 0*I, 0*I];...
         [0*I, 0*I, 0*I, 0*I,   I, 0*I];...
         [0*I, 0*I, 0*I, 0*I,  -I, 0*I]];

dJdk2 = [[0*I, 0*I, 0*I, 0*I, 0*I, 0*I];...
         [0*I, 0*I,  -I, 0*I, 0*I, 0*I];...
         [0*I, 0*I,   I, 0*I, 0*I, 0*I];...	
         [0*I, 0*I, 0*I, 0*I, 0*I, 0*I];...
         [0*I, 0*I, 0*I, 0*I, 0*I,  -I];...
         [0*I, 0*I, 0*I, 0*I, 0*I,   I]];

dJda = [[  I, 0*I, 0*I, 0*I, 0*I, 0*I];...
        [ -I, 0*I, 0*I, 0*I, 0*I, 0*I];...
        [0*I, 0*I, 0*I, 0*I, 0*I, 0*I];...	
        [0*I, 0*I, 0*I,   I, 0*I, 0*I];...
        [0*I, 0*I, 0*I,  -I, 0*I, 0*I];...
        [0*I, 0*I, 0*I, 0*I, 0*I, 0*I]];

dJdd = [[0*I,  -I, 0*I, 0*I, 0*I, 0*I];...
        [0*I,   I, 0*I, 0*I, 0*I, 0*I];...
        [0*I, 0*I, 0*I, 0*I, 0*I, 0*I];...	
        [0*I, 0*I, 0*I, 0*I,  -I, 0*I];...
        [0*I, 0*I, 0*I, 0*I,   I, 0*I];...
        [0*I, 0*I, 0*I, 0*I, 0*I, 0*I]];

dJdr = [[0*I,  -I, 0*I, 0*I, 0*I, 0*I];...
        [0*I,   I, 0*I, 0*I, 0*I, 0*I];...
        [0*I, 0*I, 0*I, 0*I, 0*I, 0*I];...	
        [0*I, 0*I, 0*I, 0*I,  -I, 0*I];...
        [0*I, 0*I, 0*I, 0*I,   I, 0*I];...
        [0*I, 0*I, 0*I, 0*I, 0*I, 0*I]];


D = mfactor(FJ,-[dJdk1*M,...
	         dJdk2*M,...
		 dJda*M,...
		 dJdd*M,...
                 dJdr*M]);

D2 = mfactor(FJ,[-2*dJdk1*D,...
	         -2*dJdk2*D,...
		 -2*dJda*D,...
		 -2*dJdd*D,...
                 -2*dJdr*D]);

