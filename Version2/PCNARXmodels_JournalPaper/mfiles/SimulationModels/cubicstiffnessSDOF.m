function xdot = cubicstiffnessSDOF(t,x,ut,u,Sys)

% cubicstiffnessTDOF simulates a 2-DOF system with a cubic stiffness (see [1])
%
% XDOT = cubicstiffnessTDOF(T,X,UT,USYS)
%
% T: time input
% X: state vector ([x1 x1dot x2 x2dot z]^T)
% UT: excitation input time index
% U: excitation input
% SYS: 2-DOF system parameters (structure)
%       SYS.m: system mass (kg)
%       SYS.k: system linear stiffness (N*s/m)
%       SYS.knl: system non-linear (cubic) stiffness (N*s/m)
%
% XDOT: State vector derivative ([xdot xddot]^T)
%
% EXAMPLE
% =========================================================================
% Sys.m = 1;       
% Sys.k = 1;  
% Sys.c = 0.001;  
% Sys.knl = 0.5;
% Tspan = 0:0.5:500;
% IC = zeros(4,1);
% ut = Tspan;
% u = 1000*sin(ut);
% [Tsim,Ysim] = ode45(@(t,x) cubicstiffnessTDOF(t,x,ut,u,Sys),Tspan,IC); 

%   Copyright 1980-2012, The MinaS Inc.
%   $ Version: 1.0 $ $Date: 15/11/2012 $

%         SDOF Simulated System
% /|   c     __________
% /|   __   |          |
% /|---__|--|          |
% /|        |          |
% /|  knl   |    m1    |---> u
% /|-/\_/\--|          |
% /|        |          |
% /|-/\/\/\-|          |
% /|   k    |__________|
% /|             |---------> x
%
%    ..    .           
%  m*x + c*x + k*x + knl*x^3 = u        

% System parameters
if nargin<5 || isempty(Sys)
    m = 1;      k = 1;      c = 0.001;      knl = 0.5;
else
    m = Sys.m;	k = Sys.k;  c = Sys.c;      knl = Sys.knl;
end

if length(ut)>1
    % Input forces
    u = interp1(ut,u,t); % Interpolate the data set (ut,u) at time t
end

% State vector                               
xdot(1,1) = x(2);                            
xdot(2,1) = (- k*x(1) - c*x(2) - knl*(x(1)^3) + u)/m;