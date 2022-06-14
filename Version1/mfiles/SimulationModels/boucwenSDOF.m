function xdot = boucwenSDOF(t,x,ut,u,BW,Sys)

% boucwenSDOF simulates a SDOF system with a hysteretic restoring force 
%   behaviour described by the Bouc-Wen model 
%
% XDOT = boucwenSDOF(T,X,UT,U,BW,SYS)
%
% T: time input
% X: state vector ([x(t) xdot(t) z(t)]^T)
% UT: excitation input time index 
% U: excitation input
% BW: Bouc-Wen model parameters (structure)
%       BW.n: parameter that controls hysteresis loop sharpness (n > 0)
%       BW.gamma: parameter that controls hysteresis loop shape (gamma in [-B,B])
%       BW.A: parameter that controls hysteresis loop scale (A = 1)
%       BW.B: parameter that controls hysteresis loop shape (B > 0)
% SYS: SDOF system parameters (structure)
%       SYS.m: system's mass (kg)
%       SYS.c: system's linear viscous damping coefficient (kg)
%       SYS.ki: system's pre-yield elastic stiffness (N/m)
%       SYS.kf: system's post-yield elastic stiffness (N/m)
%
% XDOT: State vector derivative ([xdot(t) xddot(t) zdot(t)]^T)
%
% EXAMPLE
% =========================================================================
% BW.n = 2;    BW.gamma=0.5;    BW.A=1;         BW.B = 1;
% Sys.m = 1;   Sys.c = 0.5;     Sys.ki = 600;   Sys.kf = 100;
% Tspan = 0:0.01:100;
% IC = zeros(3,1);
% ut = Tspan;
% u = 1000*sin(ut);
% [Tsim,Ysim] = ode45(@(t,x) boucwenSDOF(t,x,ut,u,BW,Sys),Tspan,IC); 
%
% [1] Song J. and Der Kiureghian A. (2006). "Generalized Bouc–Wen model for 
%       highly asymmetric hysteresis", Journal of Engineering Mechanics,
%       Vol 132, No. 6, pp. 610-618.

%   Copyright 1980-2012, The MinaS Inc.
%   $ Version: 1.0 $ $Date: 05/09/2012 $

% System parameters
if nargin<6 || isempty(Sys)
    m = 1;      c = 0.5;    ki = 600;       kf = 100;
else
    m = Sys.m;  c = Sys.c;  ki = Sys.ki;    kf = Sys.kf;
end
% Bouc-Wen Hysteresis model parameters
if nargin<5 || isempty(BW)
    n = 2;      gamma = 0.5;        A = 1;      B = 1;
else
    n = BW.n;   gamma = BW.gamma;   A = BW.A;   B = BW.B;
end

% Input force
if length(ut)>1
    u = interp1(ut,u,t); % Interpolate the data set (ut,u) at time t
else
    u = u(1);
end


% State variable
xdot(1,1) = x(2);
xdot(2,1) = -(c/m)*x(2) - (kf/m)*x(1) - ((ki-kf)/m)*x(3) + u/m;
xdot(3,1) = x(2)*( A - (B*sign(x(3)*x(2))+ gamma)*(abs(x(3))^n));