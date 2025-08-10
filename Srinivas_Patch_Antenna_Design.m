
%-------------------------------------------------------
clearvars
clc
close all
%-------------------------------------------------------
G = 10^9;
M = 10^6;
c = 10^-2;
m = 10^-3;
u = 10^-6;
n = 10^-9;
p = 10^-12;
f = 10^-15;
%-------------------------------------------------------
f0 = 16*G;
substrate_Thickness = 0.508*m;
er = 2.2; % Duroid 5880
er_t = 2.1; % Teflon 
uc = 2.99792458 *10^8; % m/s
n_0 = 376.7303; % Ohms
Z0 = 50;
%-------------------------------------------------------
lfs = uc/f0; % Freespace wavelength
ko = 2*pi / lfs; %
hmax = (0.3/(2*pi*sqrt(er)))*lfs; % maximum thickness
%-------------------------------------------------------
Print_Real_Unit('lfs',lfs,'m')
Print_Real_Unit('ko',ko,'rad/m')
Print_Real_Unit('hmax',hmax,'m')
Print_Break
%-------------------------------------------------------
% Design 01
%-------------------------------------------------------
%h = 1.5748*m; %design param
h = 0.508*m;
%-------------------------------------------------------
lambda_g = lfs/sqrt(er);
Le = (1/2)*lambda_g;
W_prime = 1.5*Le;
%-------------------------------------------------------
Er_eff = [(1/2)*(er+1)]+[((1/2)*(er-1))]*[(1+12*(h/W_prime))]^-(1/2);
delta_L = 0.412 * h * ((Er_eff + 0.3) / (Er_eff - 0.258)) * (((W_prime / h) + 0.264) / ((W_prime / h) + 0.8));
Lp = Le - 2*delta_L;
W = 1.5*Lp;
%-------------------------------------------------------
Print_Real_Unit('Er_eff',Er_eff,'')
Print_Real_Unit('delta_L',delta_L,'m')
Print_Real_Unit('Lp',Lp,'m')
Print_Real_Unit('W',W,'m')
%-------------------------------------------------------
Print_Break
Print_Real_Unit('lambda_g',lambda_g,'m')
Print_Real_Unit('Le',Le,'m')
Print_Real_Unit('W_prime',W_prime,'m')
%-------------------------------------------------------
[Gedge,G1,G12,B1,~,~]=...
    EE457_Microstrip_Patch_Conductance(Le,W,h,ko,lfs,n_0);
Rin = Z0;
Redge = 1/Gedge;
xf = (Le/pi)*acos(sqrt(Rin/Redge));
x0 = (1/2)*Le - xf;
Gedge = 2*(G1+G12);
%-------------------------------------------------------
Print_Break
Print_Real_Unit('G1',G1,'S')
Print_Real_Unit('G12',G12,'S')
Print_Real_Unit('B1',B1,'S')
Print_Real_Unit('Gedge',Gedge,'S')
Print_Real_Unit('Redge',Redge,'Ohm') % Resistance at the edge of the patch
Print_Real_Unit('xf',xf,'m')
%%
lambda0 = lfs;
f0 = 16*G;
BWf = 3.7771 * ((er - 1) / (er^2)) * (h / lambda0) * ( W/Lp) * f0;
BWp = (BWf*100)/f0;
Print_Real_Unit('BWf',BWf,'Hz')
Print_Real_Unit('BWp',BWp,'%')