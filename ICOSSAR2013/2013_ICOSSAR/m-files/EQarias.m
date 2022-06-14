function [Ia,t5,t45,t95,Iat,D5_45,D5_95] = EQarias(y,Ts)

% EQarias computes the Arias intensity, tmid, D5_45, D5_95
% criteria of a given earthquake accelererogram
N = length(y);
T = (0:N-1)/(N-1);

% Arias intensity
Ia = (pi/(2*9.81))*sum(y.^2)/N;

% Arias intensity as a function of time
Iat = (pi/(2*9.81))*cumsum(y.^2)/N;

% Normalized Arias intensity as a function of time
Iaindx = Iat./Ia;

% t5 
indx = find(Iaindx>0.05);
t5 = T(indx(1));

% t45 (tmid)
indx = find(Iaindx>0.45);
t45 = T(indx(1));

% t95
indx = find(Iaindx>0.95);
t95 = T(indx(1));

% D5_45
D5_45 = (t45-t5);

% D5_95
D5_95 = (t95-t5);

if nargout == 0
    figure(1)
    plot(T,Iaindx),hold on
    plot([t5 t5],[0 1],'--r')
    plot([t45 t45],[0 1],'--r')
    plot([t95 t95],[0 1],'--r')
    xlabel('Sample')
    ylabel('Normalized Arias intensity')
end