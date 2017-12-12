clc
clear all
close all

n234=10.5; n230=9.19E-6;
% u238=2700; u234=3078;
% dissolved thorium-230 at 313m (d0a), 524m (d0b) and 1918m (d0c). the data
% comes from Roy-Barman et al., 2002.the file stored at
% C:\Users\xuele_000\Documents\Previous Flux Data\MedFlux\MedFlux Th
% d0a=0.13; d0b=0.16; d0c=0.17;

% dissolved thorium-234 at 313m (d4a), 524m (d4b) and 1918m (d4c). the data
% Kirk's measurements. the file stored at
% C:\Users\xuele_000\Documents\Previous Flux Data\MedFlux\MedFlux 2005
% d4a=2525; d4b=2500; d4c=2520;
%
% Dab Dac Dbc represent the distance between two traps. 'a' means trap at
% 313m; 'b' means trap at 524m; 'c' means the traps at 1918m for example,
% Dab mean distance between trap at 313m and at 524m.
%
Dab=-211; Dac=-1605; Dbc=-1394;

%====================Reading data in. Assigning the file name to 'filenmae'
%facilitating reading in and writing out.

load Steady_Th234_Th230.mat
II = ones(length(Th0),1);
ji = length(II);

for n=1:5000
    
    x0 = [a*II,d*II,r*II,k1*II,k2*II,d234,Th234,...
          TH234,d230,Th230,TH230,P_s,P_l];
    
    x0 = log(x0);
    
    a    = x0(1:ji);
    d    = x0(ji+1:2*ji);
    r    = x0(2*ji+1:3*ji);
    k1   = x0(3*ji+1:4*ji);
    k2   = x0(4*ji+1:5*ji);
    d234 = x0(5*ji+1:6*ji);
    Th4  = x0(6*ji+1:7*ji);
    TH4  = x0(7*ji+1:8*ji);
    d230 = x0(8*ji+1:9*ji);
    Th0  = x0(9*ji+1:10*ji);
    TH0  = x0(10*ji+1:11*ji);
    P_s  = x0(11*ji+1:12*ji);
    P_l  = x0(12*ji+1:13*ji);
    
    var = 0.05*x0;
    Cdd = log(1+var.^2); % this formula is prescribed y Micheal & Lam
                         % (2009) on page 9 of 23, where var is relative error.
    
    Cdd = d0(Cdd);
    
    f = [[(U238-exp(d234))*n234-exp(d234)*exp(a)+exp(Th4)*(exp(d)+exp(r))];...
         [exp(d234)*exp(a)+exp(TH4)*exp(k2)-exp(Th4)*(exp(k1)+exp(r)+exp(d)+n234)];...
         [exp(Th4)*exp(k1)-exp(TH4)*(exp(k2)+n234)-PFD*exp(TH4)];...
         [(U234-exp(d230))*n230-exp(d230)*exp(a)+exp(Th0)*(exp(d)+exp(r))];...
         [exp(d230)*exp(a)+exp(TH0)*exp(k2)-exp(Th0)*(exp(k1)+exp(r)+exp(d)+n230)];...
         [exp(Th0)*exp(k1)-exp(TH0)*(exp(k2)+n230)-PFD*exp(TH0)];...
         [exp(P_l)*exp(k2)-exp(P_s)*(exp(k1)+exp(r))];...
         [exp(P_s)*exp(k1)-exp(P_l)*exp(k2)-PFD*exp(P_l)]];
    
    J11 = -d0(exp(d234)*exp(a));
    J12 = d0(exp(Th4)*exp(d));
    J13 = d0(exp(Th4)*exp(r));
    I16 = -d0(exp(d234)*(n234+exp(a)));
    I17 = d0(exp(Th4)*(exp(d)+exp(r)));
    
    J21 = d0(exp(d234)*exp(a));
    J22 = -d0(exp(Th4)*exp(d));
    J23 = -d0(exp(Th4)*exp(r));
    J24 = -d0(exp(Th4)*exp(k1));
    J25 = d0(exp(TH4)*exp(k2));
    J26 = d0(exp(d234)*exp(a));
    J27 = -d0(exp(Th4)*(exp(k1)+exp(r)+exp(d)+n234));
    J28 = d0(exp(TH4)*exp(k2));
    
    J34 = d0(exp(Th4)*exp(k1));
    J35 = -d0(exp(TH4)*exp(k2));
    J37 = d0(exp(Th4)*exp(k1));
    J38 = -d0(exp(TH4)*(exp(k2)+n234)+PFD*exp(TH4));
    
    J41 = -d0(exp(d230)*exp(a));
    J42 = d0(exp(Th0)*exp(d));
    J43 = d0(exp(Th0)*exp(r));
    J49 = -d0(exp(d230)*(n230+exp(a)));
    J410 = d0(exp(Th0)*(exp(d)+exp(r)));
    
    J51 =  d0(exp(d230)*exp(a));
    J52 = -d0(exp(Th0)*exp(d));
    J53 = -d0(exp(Th0)*exp(r));
    J54 = -d0(exp(Th0)*exp(k1));
    J55 =  d0(exp(TH0)*exp(k2));
    J59 = d0(exp(d230)*exp(a));
    J510 = -d0(exp(Th0)*(exp(k1)+exp(r)+exp(d)+n234));
    J511 = d0(exp(TH0)*exp(k2));
    
    J64 =  d0(exp(Th0)*exp(k1));
    J65 = -d0(exp(TH0)*exp(k2));
    J610 = d0(exp(Th0)*exp(k1));
    J611 = -d0(exp(TH4)*(exp(k2)+n234)+PFD*exp(TH4));
    
    % dP_sdt = [exp(P_l)*exp(k2)-exp(P_s)*(exp(k1)+exp(r))];
    J73 = -d0(exp(P_s)*exp(r));
    J74 = -d0(exp(P_s)*exp(k1));
    J75 =  d0(exp(P_l)*exp(k2));
    J712 = -d0(exp(P_s)*(exp(k1)+exp(r)));
    J713 =  d0(exp(P_l)*(exp(k2)+n234));
    
    % dP_ldt = [exp(P_s)*exp(k1)-exp(P_l)*exp(k2)-PFD*exp(P_l)]
    J84 =   d0(exp(P_s)*exp(k1));
    J85 =  -d0(exp(P_l)*exp(k2));
    J812 =  d0(exp(P_s)*exp(k1));
    J813 = -d0(exp(P_l)*(exp(k2)+n234)+PFD*exp(TH4));
    
    F  = [[J11, J12, J13, 0*I, 0*I, J16, J17, 0*I,  0*I,  0*I,  0*I,  0*I,  0*I];...
          [J21, J22, J23, J24, J25, J26, J27, J28,  0*I,  0*I,  0*I,  0*I,  0*I];...
          [0*I, 0*I, 0*I, J34, J35, 0*I, J37, J38,  0*I,  0*I,  0*I,  0*I,  0*I];...
          [J41, J42, J43, 0*I, 0*I, 0*I, 0*I, 0*I,  J49, J410,  0*I,  0*I,  0*I];...
          [J51, J52, J53, J54, J55, 0*I, 0*I, J59, J510, J511,  0*I,  0*I,  0*I];...
          [0*I, 0*I, 0*I, J64, J65, 0*I, 0*I, 0*I,  0*I, J610, J611,  0*I,  0*I];...
          [0*I, 0*I, J73, J74, 0*I, 0*I, 0*I, 0*I,  0*I,  0*I,  0*I, J712, J713];...
          [0*I, 0*I, 0*I, J84, J85, 0*I, 0*I, 0*I,  0*I,  0*I,  0*I, J812, J813]];
        
    xi=x0;
    keyboard
    %======================================total inversion, based
    %on References: Tarantola 1982a, 1982b.
    x1 = xi'+C0*F0'*inv(F0*C0*F0')*(F0*(x0-xi)'-f0');
        
    cutoff=(abs(x1-x0')<abs(x0')*0.01);
    t=length(cutoff);
    
    if all(cutoff)
        
        disp(n);
        
        break;
        
    else
        
        x0=x1';
        
    end
    
    if (rem(n,5)==0);
        
        disp({'the iteration time is:' n});
        disp(cutoff);
        
    end
    
end


C=C0-C0*F0'*inv(F0*C0*F0')*F0*C0;
      

      



