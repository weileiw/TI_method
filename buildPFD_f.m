function [PFdiv,dPFDdk2] = buildPFD_f(M3d,p,grd);
% this code is used to build Particle Flux Diverngence (PFD)
%
%                      ________________________                                        
%                     /         A             /|  POP sinking          
%                 top/_______________________/ |       |       
%                    |  |                    | |       | -w   
%                    | /        V            | /       |       
%                 bot|/______________________|/        V
%                                                  
% PFD = (A*w(top)*POP(top)-A*w(bot)*POP(bot))/V;
% add an exra layer of zeros at the bottom to ensure there is no
% flux in or out of the domain when using periodic shift operators
% for finite differences and averaging
spd = 24*60*60;
[ny,nx,nz] = size(M3d);
M3D = zeros(ny,nx,nz+1);
M3D(:,:,1:end-1) = M3d;
% add the zw coordinate at the top of the extra layer
ZW3d = grd.ZW3d;
ZW3d = ZW3d(:,:,[1:end,end]);
ZW3d(:,:,end) = grd.ZW3d(:,:,end)+grd.dzt(end);
% areas of the top of the grid box
dAt = (grd.DXT3d.*grd.DYT3d).*M3d;
% volume of the grid boxes
dVt = (grd.DXT3d.*grd.DYT3d.*grd.DZT3d).*M3d;
%
n = nx*ny*(nz+1);
I0 = speye(n);
i0 = zeros(ny,nx,nz+1);
i0(:) = 1:n;
% periodic shifts OK because M3D has a layer of zeros on the bottom
iu = i0(:,:,[nz+1,1:nz]); %use a periodic upward shift
ib = i0(:,:,[2:nz+1,1]); % use a periodic downward shift
IU = I0(iu,:);
IB = I0(ib,:);
% keep only wet boxes
iocn = find(M3D(:));
I0 = I0(iocn,:); I0 = I0(:,iocn);
IU = IU(:,iocn); IU = IU(iocn,:);
IB = IB(:,iocn); IB = IB(iocn,:);
% (averages POP onto the top of the grid boxes)
AVG = d0((I0+IU)*squeeze(M3D(iocn)))\(I0+IU);
% (compute the divergence in the center of the grid boxes)
DIV = d0(dVt(iocn))\(I0-IB)*d0(dAt(iocn));
% (compute the flux at the top of the grid cells)
% mimics a Martin curve flux attenuation profile
%(see Kriest and Oschelies 2008 in Biogeosciences)
%r = p.kappa_p;

r = p.k2;
b = p.bf;
a = r/b;

% particle sinking velocity at the top of the grid cells.
MSK = M3D.*M3D(:,:,[nz+1,1:nz]);
M = MSK.*ZW3d;
w = -a*M;
%w = -100;
dadb = -r/(b^2);
dadk2 = 1/b;
dwdb = -dadb*M;
dwdk2 = -dadk2*M;
d2adb2 = 2*r/(b^3);
d2wdb2 = -d2adb2*M;
FLUX = d0(w(iocn))*AVG;
dFLUXdb = d0(dwdb(iocn))*AVG;
dFLUXdk2 = d0(dwdk2(iocn))*AVG;
d2FLUXdb2 = d0(d2wdb2(iocn))*AVG;
%FLUX = d0(w(iocn))*IU;
%dFLUXdb = d0(dwdb(iocn))*IU;
%dFLUXdk2 = d0(dwdk2(iocn))*IU;
%d2FLUXdb2 = d0(d2wdb2(iocn))*IU;
% particle flux divergence operator
PFdiv = DIV*FLUX;
dPFDdb = DIV*dFLUXdb;
dPFDdk2 = DIV*dFLUXdk2;
d2PFDdb2 = DIV*d2FLUXdb2;

end 
