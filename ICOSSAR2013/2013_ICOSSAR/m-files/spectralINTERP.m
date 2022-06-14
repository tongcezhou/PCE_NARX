function [th,tmid] = spectralINTERP(Y,Ts,options)

N = length(Y);
Td = (0:N-1);
Tc = Ts*(0:N-1);
% +Zero crossings
ZC = [0 ;cumsum(Y(2:end)>0 & Y(1:end-1)<0)];

% Arias intensity of the real EQ accelerogram
[Ia,~,~,~,Iat] = EQarias(Y,Ts); 
% Normalized Arias intensity as a function of time
Iaindx = Iat./Ia;

% t1 
indx = find(Iaindx>0.01);
t1 = Td(indx(1));
% t99
indx = find(Iaindx>0.99);
t99 = Td(indx(1));

% X points on the ZC curve between 1 and 99% of the Arias intensity function
Tint = round(linspace(t1,t99,options.points));
% Fit second order polynomial
p = polyfit(Tc(Tint),ZC(Tint)',2);

% tmid
indx = find(Iaindx>0.45);
tmid = Td(indx(1));
% omega_mid
th(1) = Tc(tmid).*p(1)+ p(2);
% omega_dot
th(2) = p(1);