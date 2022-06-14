function [Yp,E] = narxpred(Y,X,n,regstr,theta)

%   NARXPRED Computes the one-step-ahead predictions of a NARX model for a
%   given signal:
%   ^
%   y[t|t-1] =  - \sum_{i=1}^{Nr} theta_i * f(y[t-j], x[t-k])
%
%	[Yp,E] = NARXPRED(Y,X,T0,REGSTR,THETA)
%
%   Y : Output data (T x 1)
%   X : Excitation input data (T x 1)
%   T0 : Initial time (max lag in the rergessors + 1)
%   REGSTR : Cell array consisting of the regressors in string format e.g.
%               '(Y(t-2,:).^3)*(X(t,:))'
%   THETA : NARX coefficient parameter vector (Nr x 1)
%
%   Yp : One-step-ahead predictions (T x 1)
%   E : estimated residual series (T x 1)

%   Copyright 1980-2013, The MinaS Inc.
%   $ Version: 1.01 $ $Date: 19/04/2013 $

% Number of samples
T = size(X,1);
% Construction of the regressor matrix
str = '';
for j = 1:length(regstr)
    str = strcat(str,'phi(',num2str(j),')=',regstr{j},';');
end

% Calculation of the one-step ahead predictions
Yp = zeros(T,1);

for t = n+1:T
    eval(str);
    Yp(t) = theta'*phi';
end
E = Y - Yp;