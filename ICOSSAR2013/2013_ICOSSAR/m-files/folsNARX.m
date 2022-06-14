function [nr,REGstr,PHIind,res,criteria] = folsNARX(Y,X,orders,regstr,degree,q,maxterms,Nr,thr,options)

%   folsNARX Forward Orthogonal Least Squares NARX model structure selection
%
%   y[t] + \sum_{i=1}^{Nr} theta_i * f(y[t-j], x[t-k]) = e[t]
%
%	[Nr,REGSTR,PHIind,RES,CRITERIA] = folsNARX(Y,X,ORDERS,REGSTR,P,Q,MAXTERMS,Nr,THR,OPTIONS)
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
%   Nr : Number of regressors to be selected for the final model
%   THR : The threshold, which has to be a decimal number in the interval
%           [0,1], indicates the maximum relative error allowed between the
%           approximated output yr and the desired output y, when the
%           number of regressors used is r.
%   OPTIONS : estimation method options (structure)
%             criterion: 'err' or 'bic'
%             t0 : number of residuals excluded from criteria computation
%
%   Nr : Number of regressors finally selected
%   REGSTR : Cell array consisting of the finally selected regressors
%   PHIind : Index of the finally selected regressors in the initial regressor matrix 
%   RES : estimated residual series (T x 1)
%   CRITERIA :
%   	SPP : Samples Per Parameter
%   	CONNUM : Regressor matrix condition number
%       RSS : Residual Sum of Squares
%       RSS_SSS : Normalized Residual Sum of Squares
%       BIC : Bayesian Information Criterion

%   * If REGSTR is empty the regressor terms are formed by polynomial basis
%       functions of maximum total degree equal to P, maximum quasi-norm Q,
%       and maximum number of terms in each regressor equal to MAXTERMS.

%   Copyright 1980-2013, The MinaS Inc.
%   $ Version: 1.1 $ $Date: 18/04/2013 $

if ~(isfield(options,'t0')),                            options.t0 = 100;           end
if ~(isfield(options,'criterion')),                     options.criterion = 'rss';  end
if ~(isfield(options,'Nr_step')),                       options.Nr_step = 20;       end
if strcmp(options.criterion,'bic') && isempty(Nr),      Nr = 100;                   end

% Number of samples
T = size(Y,1);
% Number of initial residuals excluded from the PE criterion
t0 = options.t0;

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
Nrmax = length(regstr);

% Regressor matrix construction
t = n+1:T;
PHI = zeros(T-n,length(regstr));
for j = 1:length(regstr)
    PHI(1:T-n,j) = eval(regstr{j});
end

Ystack = Y(n+1:end,:);

% PHIind : Index of selected regressors
% ActReg : Active regressors
% PHIo : Finally selected regressor matrix
% THo : Coefficients for the selected regressors

% The algorithm is terminated when the added regressors entail increment
% of the RSS
if strcmp(options.criterion,'bic')
    flag = 0;
    while flag == 0
        BIC0 = inf;
        [~,PHIind,~,PHIo,THo] = fols(PHI,Ystack,Nr);
        
        for nr = 1:Nr
            % Residuals computation
            E = [zeros(n,1); Ystack - (PHIo(:,1:nr)*THo(1:nr))];
            est_param = nr;
            % BIC computation        
            BIC = log(var(E(t0+1:T))) + est_param*log(T-t0)/(T-t0);              
            if BIC < BIC0
                criteria.bic_history(nr) = BIC;
                BIC0 = BIC;
            else
                flag = 1;
                break
            end
        end
        if (flag ==0) && (Nr + options.Nr_step < Nrmax);
            Nr = Nr + options.Nr_step;
        else
            Nr = Nrmax;
        end
    end
    
    % Phi may be truncated to nr-1 terms (inclusion of the nr-th term leads
    %   to BIC criterion increment.
    nr = nr-1;
    PHIind = PHIind(1:nr);
    PHIo = PHIo(:,1:nr);
    THo = THo(1:nr);
    REGstr = cell(1,nr);
    for ii = 1:nr
        REGstr{1,ii} = regstr{PHIind(ii)};
    end
    
% ERR criteria is utilized for the termination of the algorithm
else
    [~,PHIind,~,PHIo,THo] = fols(PHI,Ystack,Nr,thr);
    
    nr = length(PHIind);
    REGstr = cell(1,nr);
    for ii = 1:nr
        REGstr{1,ii} = regstr{PHIind(ii)};
    end
end

% Residuals computation
res = Ystack - (PHIo*THo);
criteria.rss_sss = 100*norm(res)^2/norm(Ystack)^2;
if criteria.rss_sss < inf
    est_param = nr;
    criteria.spp = (T-n)/est_param;
    criteria.connum = cond(PHIo);
    criteria.bic = log(var(res)) + est_param*log(T-n)/(T-n);
else
    est_param = nr;
    criteria.spp = (T-n)/est_param;
    criteria.connum = cond(PHIo);
    criteria.bic = Inf;
end