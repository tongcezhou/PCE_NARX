function [thetaij,res,criteria] = Runpad_narx(Y,X,Q,orders,indx,customreg,th0,options)

%   PC-NARX Computes parameter estimates for an Polynomial Chaos
%   AutoRegressive model with eXogenous input:
%
%   y[t] + \sum_{i=1}^{Nr} a_i(xi) * f(y[t-j],x[t-k]) =  e[t]
%
%	[THETAij,RES,CRITERIA] = PCNARX(Y,X,Q,ORDERS,INDX,CUSTOMREG,TH0,OPTIONS)
%
%   Y : Output data (T x K)
%   X : Excitation input data (T x K)
%   Q : Input variables data (M x K)
%   ORDERS = [Na Nd Nb] : The AR/X orders of the PC-ARX input-output model
%   INDEX : basis function indices (cell)
%       INDX{1} : d(j) basis function multi-indices vector
%   OPTIONS : estimation method options (structure)
%       focus : 'prediction'/'simulation'
%       basis : 'hermi'/'lague'/'cheby'/'legen'/'jacob' (see PCbasis function)
%       t0 : number of residuals excluded from criteria computation
%
%   THETAij : projection coefficient matrix (Nr*P x 1)
%   E : estimated residual series (T x K)
%   CRITERIA :
%   	SPP : Samples Per Parameter
%       RSS : Residual Sum of Squares
%       RSS_SSS : Normalized Residual Sum of Squares
%       BIC : Bayesian Information Criterion

%   Copyright 1980-2013, The MinaS Inc.
%   $ Version: 1.11 $ $Date: 18/04/2013 $

global T M K n p B Nr basis maxsize

% Polynomial basis option
if ~(isfield(options,'basis')),     options.basis = 'legen';        end
% Estimation method option
if ~(isfield(options,'focus')),     options.focus = 'prediction';  	end
% Initial values for NL optimization (focus simulation) based on PE method
if ~(isfield(options,'PEinit')),    options.PEinit = 'y';         	end
% Allowable matrix size (for larger regression matrix OLS is implemented by sums)
if ~(isfield(options,'maxsize')),   options.maxsize = 1e9;          end
maxsize = options.maxsize;

% Transformation of the input vector (mormalized to the standard distribution)
if isfield(options,'TransM')
    Q = eval(options.TransM);
end

% Number of samples
T = size(Y,1);
% Number of variables
M = size(Q,1);
% Number of experiments
K = size(Q,2);
% Number of regressors
Nr = length(customreg);
% Maximum univariate polynomial degree in the PC basis
p = max(max(indx{1}));
% Total number of basis functions
B = size(indx{1},1);
% Basis constraction: size(basis_i) = p x K
basis = cell(M,1);
for m = 1:M
    basis{m} = PCbasis(Q(m,:),p,options);
end
stop

% PC-ARX case
    % AR order
    na = orders(1);
    % X order
    nb = orders(3);
    % Initial time instant
    n = max(na,nb);
    if strcmp(options.focus,'prediction')
        [thetaij,res,criteria] = pcnarxPE(Y,X,indx,customreg);
    else
        if isempty(th0) && ( ~(isfield(options,'PEinit')) || strcmp(options.PEinit,'n') )
            th0 = zeros(Nr*B,1);
        else
            th0 = narxPE(Y,X,indx,customreg);
        end
        [thetaij,res,criteria] = pcnarxSIM(th0,Y,X,indx,customreg,options);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Nonlinear PE Refinement
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [THETA,res,criteria] = narxPE(Y,X,indx,customreg)

global T M K Nr n B basis maxsize 

% Solution using the summation method (not-recommended).
% !!! Use this method only when the whole regression matrix
%     can't be constructed due to memory problems
if K*(T-n)*Nr*B > maxsize
    
    disp(['The regressor matrix has more than ',num2str(maxsize),' elements.',...
        ' Estimation method is based on sums.'])
    % Maximum number of columns of the regression matrix
    Tdata = floor(maxsize/(K*Nr*B));
    % Number of iterations
    Titer = ceil((T-n)/Tdata);
    
    PHI_PHI = zeros(Nr*B);
    PHI_Y = zeros(Nr*B,1);
    % Regression matrix formulation (AR basis part)
    phiB = ones(B,K);
    for i = 1:B
        for k = 1:K
            for m = 1:M
                phiB(i,k) = phiB(i,k)*basis{m}(indx{1}(i,m)+1,k);
            end
        end
    end
    phiY = kron(ones(1,Nr),phiB');
    clear phiB
    BregY = kron(ones(Tdata,1),phiY);
   
    % Iterations 1:Titer-1
    for i = 1:Titer-1
        Ystack = reshape(Y(n+1+(i-1)*Tdata:n+i*Tdata,:)',K*Tdata,1);
        % Y regressor part
        Yreg = zeros(K*Tdata,n);
        t = n+1+(i-1)*Tdata:n+i*Tdata;
        for j = 1:Nr
            Yreg(1:K*Tdata,j) = reshape(eval(customreg{j})',K*Tdata,1);
        end
        Yreg = kron(Yreg,ones(1,B));
        PHI = Yreg.*BregY;
        
        PHI_PHI = PHI_PHI + PHI.'*PHI;
        PHI_Y = PHI_Y + PHI.'*Ystack;
    end

    BregY = kron(ones(T-(Titer-1)*Tdata-n,1),phiY);
    % Final iteration
    Ystack = reshape(Y(n+1+(Titer-1)*Tdata:T,:)',K*(T-n-(Titer-1)*Tdata),1);
    % Y regressor part
    Yreg = zeros(K*(T-n-(Titer-1)*Tdata),n);
    t = n+1+(Titer-1)*Tdata:T;
    for j = 1:Nr
        Yreg(1:K*(T-n-(Titer-1)*Tdata),j) = reshape(eval(customreg{j})',K*(T-n-(Titer-1)*Tdata),1);
    end
    Yreg = kron(Yreg,ones(1,B));
    PHI = Yreg.*BregY;
    PHI_PHI = PHI_PHI + PHI.'*PHI;
    PHI_Y = PHI_Y + PHI.'*Ystack;

    % Parameters matrix computation
    THETA = (PHI_PHI./(T-n))\(PHI_Y./(T-n));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %           Residuals computation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    res = zeros(K*(T-n),1);
    SSS = 0;
    % Iterations 1:Titer-1
    BregY = kron(ones(Tdata,1),phiY);
    for i = 1:Titer-1
        Ystack = reshape(Y(n+1+(i-1)*Tdata:n+i*Tdata,:)',K*Tdata,1);
        SSS = SSS + norm(Ystack)^2;
        % Y regressor part
        Yreg = zeros(K*Tdata,n);
        t = n+1+(i-1)*Tdata:n+i*Tdata;
        for j = 1:Nr
            Yreg(1:K*Tdata,j) = reshape(eval(customreg{j})',K*Tdata,1);
        end
        Yreg = kron(Yreg,ones(1,B));
        PHI = Yreg.*BregY;
        res((i-1)*Tdata*K+1:i*Tdata*K,1) = Ystack - PHI*THETA;
    end

    % Final iteration
    BregY = kron(ones(T-(Titer-1)*Tdata-n,1),phiY);
    Ystack = reshape(Y(n+1+(Titer-1)*Tdata:T,:)',K*(T-n-(Titer-1)*Tdata),1);
    SSS = SSS + norm(Ystack)^2;
    % Y regressor part
    Yreg = zeros(K*(T-n-(Titer-1)*Tdata),n);
    t = n+1+(Titer-1)*Tdata:T;
    for j = 1:Nr
        Yreg(1:K*(T-n-(Titer-1)*Tdata),j) = reshape(eval(customreg{j})',K*(T-n-(Titer-1)*Tdata),1);
    end
    Yreg = kron(Yreg,ones(1,B));
    PHI = Yreg.*BregY;
    res((Titer-1)*Tdata*K+1:(T-n)*K,1) = Ystack - PHI*THETA;
    clear PHI

    % Criteria Computation
    %%%%%%%%%%%%%%%%%%%%%%
    est_param = Nr*B;
    criteria.connum = sqrt(cond(PHI_PHI));
    criteria.spp = (T-n)*K/est_param;
    criteria.rss = norm(res)^2;
    criteria.rss_sss = 100*criteria.rss/SSS;
    criteria.bic = log(var(res)) + est_param*log(T-n)/(T-n);  
    % QR or SVD method
    % !!! recommended
else
    Ystack = reshape(Y(n+1:T,:)',K*(T-n),1);
    % Number of regressors
    t = n+1:T;
    % Y regressor part
    Yreg = zeros(K*(T-n),Nr);
    for j = 1:Nr
        Yreg(1:(T-n)*K,j) = reshape(eval(customreg{j})',(T-n)*K,1);
    end
    Yreg = kron(Yreg,ones(1,B));

    % Regression matrix formulation (basis part)
    Breg = ones(B,K);
    for i = 1:B
        for k = 1:K
            for m = 1:M
                Breg(i,k) = Breg(i,k)*basis{m}(indx{1}(i,m)+1,k);
            end
        end
    end
    BregY = kron(ones(1,Nr),Breg');
    BregY = kron(ones(T-n,1),BregY);
    PHI = Yreg.*BregY;

    % Solution using the QR method
    THETA = PHI\Ystack;
    % Residuals computation
    res = Ystack - PHI*THETA;

    % Criteria Computation
    %%%%%%%%%%%%%%%%%%%%%%
    est_param = Nr*B;
    criteria.spp = (T-n)*K/est_param;
    criteria.connum = cond(PHI);
    criteria.rss = norm(res)^2;
    criteria.rss_sss = 100*criteria.rss/norm(Ystack)^2;
    criteria.bic = log(var(res)) + est_param*log(T-n)/(T-n);  
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%			Nonlinear refinement based on simulation error
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [TH,res,criteria] = narxSIM(TH0,Y,X,indx,customreg,options)

global T K Nr n M B basis maxsize    
    
% Nonlinear optimization algorithm options
if ~isfield(options,'nlopts')
    if (T-n)*K > maxsize
        options.nlopts = optimset('LargeScale','off','Display','iter',...
            'TolFun',1e-6,'TolX',1e-9,'MaxFunEval',200000,'MaxIter',200000);
    else
        options.nlopts = optimset('Algorithm','Levenberg-Marquardt','Display','iter',...
            'TolFun',1e-6,'TolX',1e-9,'MaxFunEval',200000,'MaxIter',200000);
    end
end


phiB = ones(B,K);
for i = 1:B
    for k = 1:K
        for m = 1:M
            phiB(i,k) = phiB(i,k)*basis{m}(indx{1}(i,m)+1,k);
        end
    end
end

% Nonlinear optimization
if T*K > maxsize
    [TH,~,criteria.PEexitflag,criteria.PEoutput] = ...
        fminunc(@(theta0) pcnarxSIMnlopt(theta0,Y,X,phiB,customreg),TH0,options.nlopts);
else
    [TH,~,~,criteria.PEexitflag,criteria.PEoutput] = ...
        lsqnonlin(@(theta0) pcnarxSIMnlopt(theta0,Y,X,phiB,customreg),TH0,[],[],options.nlopts);
end


maxsize = inf;
res = pcnarxSIMnlopt(TH,Y,X,phiB,customreg);
res = reshape(res(1:T,:)',K*T,1);
% Criteria Computation
%%%%%%%%%%%%%%%%%%%%%%
Ystack = reshape(Y(n+1:T,:)',K*(T-n),1);
SSS = norm(Ystack(:))^2;
est_param = Nr*B;
criteria.spp = (T-n)*K/est_param;
criteria.rss = norm(res)^2;
criteria.rss_sss = 100*criteria.rss/SSS;
criteria.bic = log(var(res)) + est_param*log(T-n)/(T-n);  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%			Nonlinear refinement based on simulation error
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = pcnarxSIMnlopt(theta0,Y,X,phi,customreg)

global T K Nr B n maxsize

an = phi'*reshape(theta0,B,Nr);

Ysim = zeros(T,K);
for t = n+1:T
    for j = 1:Nr
        Ysim(t,:) = Ysim(t,:) + an(:,j)'.*eval(customreg{j});
    end
end

if T*K > maxsize
    varargout{1} = norm(Y(:)-Ysim(:))^2;
else
    varargout{1} = Y-Ysim;
end