function Y = narxsim(X,Y0,n,regstr,theta,p)

%   NARXSIM Computes the simulations of a NARX model for a given input signal:
%   ~                                     ~
%   y[t] =  - \sum_{i=1}^{Nr} theta_i * f(y[t-j], x[t-k])
%
%	Ys = NARXSIM(X,Y0,T0,REGSTR,THETA)
%
%   X : Excitation input data (T x 1)
%   Y0 : Initial output values (T0 x 1)
%   n : Initial time (max lag in the rergessors + 1)
%   REGSTR : Cell array consisting of the regressors in string format e.g.
%               '(Y(t-2,:).^3)*(X(t,:))'
%   THETA : NARX coefficient parameter vector (Nr x 1)
%
%   Ys : Simulations (T x 1)

%   Copyright 1980-2013, The MinaS Inc.
%   $ Version: 1.01 $ $Date: 19/04/2013 $

% Number of samples
T = size(X,1);

% Construction of the regressor matrix
phi = zeros(length(regstr),1);
str = '';
for j = 1:length(regstr)
    str = strcat(str,'phi(',num2str(j),')=',regstr{j},';');
end


if ~isempty(findstr(str,'theta1_0'))
    THnl = theta(length(regstr)+1:end);
    theta = theta(1:length(regstr));
    cntr = 0;
    for k = 1:length(THnl)/p
        for j = 0:p-1
            cntr = cntr+1;
            eval(['theta',num2str(k),'_',num2str(j),' = THnl(cntr);']);
        end
    end
end


% Calculation of the one-step ahead predictions
Y = zeros(T,1);
Y(1:n) = Y0;
for t = n+1:T
    eval(str);
    Y(t) = phi.'*theta;
end