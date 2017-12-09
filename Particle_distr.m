% this function intends to generate idealized data for 
% two tpye of particles in the ocean. One type is small
% particles and not sinking, the other type is large 
% particles that sink at a speed of w = 100m/d. Small particles
% reminieralized with a constant of beta_neg_1 = 1 y^-1;
% small particles also lose via particle aggregation beta_2 = 3 y^-1;
% large particles disaggregate to form small particles beta_nega_2 = 150
% y^-1;
% T(P_s) = beta_nega2*P_l - (beta_nega1+beta2)*P_s  (1);
% T(P_l) = beta2*P_s - beta_nega2*P_l - w(dP_l/dz)   (2);
function [P_l,P_s] = particle_cycl(p)
dpa        = 365;     % day per year
w          = 150*dpa; % m y^-1;
beta_nega1 = 1;       % y^-1;
beta2      = 3;       % y^-1;
beta_nega2 = 150;     % y^-1;

% depth 
z = linspace(110,4900,100);
% large particle concentration at 110m 
P_l0 = 1e-6;

lambda = -beta_nega1*beta_nega2/(w*(beta_nega1+beta2));
P_l = P_l0*exp(lambda*z);
P_s = beta_nega2*P_l./(beta_nega1+beta2);

