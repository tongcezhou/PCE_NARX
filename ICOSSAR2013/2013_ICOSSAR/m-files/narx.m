function [theta,res,criteria,regstr] = narx(Y,X,orders,regstr,degree,q,maxterms,th0,options)

%   NARX Computes parameter estimates for a Nonlinear AutoRegressive model 
%   with eXogenous input:
%
%   y[t] + \sum_{i=1}^{Nr} theta_i * f(y[t-j], x[t-k]) = e[t]
%
%	[THETA,RES,CRITERIA,REGSTR] = NARX(Y,X,ORDERS,REGSTR,P,Q,MAXTERMS,TH0,OPTIONS)
%
%   Y : Output data (T x 1)
%   X : Excitation input data (T x 1)
%   ORDERS = [Na Nd Nb] : The AR/X orders of the NARX input-output model
%   REGSTR : Cell array consisting of the regressors in string format e.g.
%               '(Y(t-2,:).^3)*(X(t,:))'
%   P : Maximum total polynomial degree for each regressor, e.g. for
%               '(Y(t-1,:).^i)*(X(t,:).^j)'  --->  i+j < P
%   Q : quasi-norms order (help combinations)
%   MAXTERMS : Maximum number of terms included in each regressor
%   TH0 : Initial parameter vector (for simulation error estimation)
%   OPTIONS : estimation method options (structure)
%       focus : 'prediction'/'simulation'
%       t0 : number of residuals excluded from criteria computation
%
%   THETA : NARX coefficient parameter vector (Nr x 1)
%   RES : estimated residual series (T x 1)
%   CRITERIA :
%   	SPP : Samples Per Parameter
%   	CONNUM : Regressor matrix condition number
%       RSS : Residual Sum of Squares
%       RSS_SSS : Normalized Residual Sum of Squares
%       BIC : Bayesian Information Criterion
%
%   * If REGSTR is empty the regressor terms are formed by polynomial basis
%       functions of maximum total degree equal to P, maximum quasi-norm Q,
%       and maximum number of terms in each regressor equal to MAXTERMS.

%   Copyright 1980-2013, The MinaS Inc.
%   $ Version: 1.1 $ $Date: 18/04/2013 $

global T n

% Estimation method option
if ~(isfield(options,'focus')),     options.focus = 'simulation';  	end
% Initial values for NL optimization (focus simulation) based on PE method
if ~(isfield(options,'PEinit')),    options.PEinit = 'y';         	end

% Number of samples
T = size(Y,1);

% NAR case
if isempty(X)
    % AR order
    na = orders(1);
    % Initial time instant
    n = na;
    % X part variables
    [nd,nb] = deal([]);
% NARX case
else
    % AR order
    na = orders(1);
    % Input delay order
    nd = orders(2);
    % X order
    nb = orders(3);
    % Initial time instant
    n = max(na,nb);
end

% Regressors 
if isempty(regstr)
    regstr = polyregNARX([na nd nb],degree,q,maxterms);
end
    
if strcmp(options.focus,'prediction') || isempty(X)
    [theta,res,criteria] = narxPE(Y,X,regstr);
else
    if isempty(th0) && ( ~(isfield(options,'PEinit')) || strcmp(options.PEinit,'n') )
        th0 = zeros(length(regstr),1);
    else
        th0 = narxPE(Y,X,regstr);
    end
    [theta,res,criteria] = narxSIM(th0,Y,X,regstr,options);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Nonlinear PE Refinement
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [THETA,res,criteria] = narxPE(Y,X,regstr)

global T n

% Regressor matrix construction
t = n+1:T;
PHI = zeros(T-n,length(regstr));
for j = 1:length(regstr)
    PHI(1:T-n,j) = eval(regstr{j});
end

% Solution using the QR method
THETA = PHI\Y(n+1:T,:);
% Residuals computation
res = Y(n+1:T,:) - PHI*THETA;

% Criteria Computation
%%%%%%%%%%%%%%%%%%%%%%
est_param = length(THETA);
criteria.spp = (T-n)/est_param;
criteria.connum = cond(PHI);
criteria.rss = norm(res)^2;
criteria.rss_sss = 100*criteria.rss/norm(Y(n+1:T,:))^2;
criteria.bic = log(var(res)) + est_param*log(T-n)/(T-n);  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%			Nonlinear refinement based on simulation error
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [THETA,res,criteria] = narxSIM(TH0,y,x,regstr,options)

global T n
    
% Nonlinear optimization algorithm options
if ~isfield(options,'nlopts')
    options.nlopts = optimset('Algorithm','Levenberg-Marquardt','Display','iter',...
        'TolFun',1e-6,'TolX',1e-9,'MaxFunEval',200000,'MaxIter',200000);
end

% Nonlinear optimization
[THETA,~,~,criteria.PEexitflag,criteria.PEoutput] = ...
    lsqnonlin(@(theta0) narxSIMnlopt(theta0,y,x,regstr),TH0,[],[],options.nlopts);

res = narxSIMnlopt(THETA,y,x,regstr);

% Criteria Computation
%%%%%%%%%%%%%%%%%%%%%%
SSS = norm(y(n+1:T,:))^2;
est_param = length(THETA);
criteria.spp = (T-n)/est_param;
criteria.rss = norm(res)^2;
criteria.rss_sss = 100*criteria.rss/SSS;
criteria.bic = log(var(res)) + est_param*log(T-n)/(T-n);  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%			Nonlinear refinement based on simulation error
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = narxSIMnlopt(theta0,Y,X,regstr)

global T n

Ysim = zeros(T,1);
for t = n+1:T
    for j = 1:length(regstr)
        Ysim(:,t) = Ysim(:,t) - theta0(j)*eval(regstr{j});
    end
end
varargout{1} = Y-Ysim;