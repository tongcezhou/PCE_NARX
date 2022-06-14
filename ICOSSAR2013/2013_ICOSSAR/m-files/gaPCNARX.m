function [INDX,thetaij,res,criteria,GAoutput] = gaPCNARX(Y,X,Q,orders,customreg,P,q,options)

%   gaPCNARX: GA-based PC functional subspace selection for a Polynomial 
%       Chaos Nonlinear AutoRegressive with eXogenous input model:
%
%   y[t] + \sum_{i=1}^{Nr} theta_i(xi) * f(y[t-j], x[t-k]) = e[t]
%
% [INDX,THETAij,RES,CRITERIA,GAoutput] = gaPCNARX(Y,X,Q,ORDERS,REGSTR,P,q,OPTIONS)
%
%   Y : Output data (T x K)
%   X : Excitation input data (T x K)
%   Q : Input variables data (M x K)
%   ORDERS = [Na Nd Nb] : The AR/X orders of the PC-NARX input-output model
%   P : Total maximum polynomial degree for the PC basis
%   q : maximum quasi-norm order
%   OPTIONS : estimation method options (structure)
%       basis : 'hermi'/'lague'/'cheby'/'legen'/'jacob' (see PCbasis function)
%       t0 : number of residuals excluded from criteria computation
%   OPTIONS.GA : Genetic algorithm options (structure)
%
%   INDX : The finally selected truncation set
%   THETAij : NARX projection coefficient vector (Nr*B x 1) where B the
%             length of the selected PC functional subspace
%   E : estimated residual series (T x K)
%   CRITERIA :
%   	SPP : Samples Per Parameter
%       RSS : Residual Sum of Squares
%       RSS_SSS : Normalized Residual Sum of Squares
%       BIC : Bayesian Information Criterion
%   OUTPUT : returns a structure OUTPUT with the following information:
%               randstate: State of the function RAND used before GA started
%               randnstate: State of the function RANDN used before GA started
%               generations: Total generations, excluding HybridFcn iterations
%               funccount: Total function evaluations
%               maxconstraint: Maximum constraint violation, (if any)
%               message: GA termination message

%   Copyright 1980-2013, The MinaS Inc.
%   $ Version: 1.01 $ $Date: 18/04/2013 $

% Global variables
global T M K n Nr B basis maxsize

% Method options
if ~(isfield(options,'warnings')),      options.warnings = 'on' ;  end
if ~(isfield(options,'basis')),         options.basis = 'legen' ;  end
if ~(isfield(options,'ortho')),         options.ortho = 'y';       end
if ~(isfield(options,'maxsize')),       options.maxsize = 1e7;     end
maxsize = options.maxsize;

if strcmp(options.warnings,'off')
    % Turning off warning for rank deficiency
    warning('off','MATLAB:rankDeficientMatrix');
    warning('off','MATLAB:nearlySingularMatrix');
end

% Genetic Algorithm options
if ~(isfield(options.GA,'PopulationSize')),      options.GA.PopulationSize = 100 ;          end
if ~(isfield(options.GA,'Generations')),         options.GA.Generations = 100;              end
if ~(isfield(options.GA,'CrossoverFraction')),   options.GA.CrossoverFraction = 0.8;        end
if ~(isfield(options.GA,'Display')),             options.GA.Display = 'iter';               end
if ~(isfield(options.GA,'PlotFcns')),            options.GA.PlotFcns = {@gaplotbestf,@gaplotbestindiv};  end
if ~(isfield(options.GA,'StallGenLimit')),       options.GA.StallGenLimit = 20;             end
if ~(isfield(options.GA,'EliteCount')),
    options.GA.EliteCount = round(0.05*options.GA.PopulationSize);                          end
if ~(isfield(options.GA,'MutationFcn')),
    options.GA.MutationFcn = {@mutationuniform,0.2};                                        end

options.GA.PopulationType = 'bitstring';
options.GA.Vectorized = 'on';
options.GA.StallTimeLimit = Inf;
GAoptions = gaoptimset(options.GA);

% Transformation of the input vector
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

% PC-AR case
if isempty(X)
    % AR order
    na = orders(1);
    % Initial time instant
    n = na;
% PC-ARX case
else
    % AR order
    na = orders(1);
    % X order
    nb = orders(3);
    % Initial time instant
    n = max(na,nb);
end

% Index of the polynomial basis
indx0{1} = combinations(M,P,q);
% Total number of basis functions
B = size(indx0{1},1);
% Basis constraction: size(basis_i) = p x K
basis = cell(M,1);
for m = 1:M
    basis{m} = PCbasis(Q(m,:),P,options);
end

% gene length
gene = size(indx0{1},1);

% Computation of the Maximum Dimensions Information Matrix
RegMat = PCARXinfmat(Y,X,indx0,customreg);

% Fitness function computation
[indiv,GAoutput.fval,GAoutput.reason,GAoutput.output,GAoutput.population,GAoutput.scores] = ...
    ga(@(gene) gaPCNARXff(gene,RegMat),gene,[],[],[],[],[],[],[],GAoptions);
INDX{1} = indx0{1}(indiv == 1,:);

[thetaij,res,criteria] = pcnarx(Y,X,Q,orders,INDX,customreg,[],options);

% Turning on warning for rank deficiency
warning('on','MATLAB:rankDeficientMatrix');
warning('on','MATLAB:nearlySingularMatrix');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              Information Matrix Computation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function RegMat = PCARXinfmat(Y,X,indx0,customreg)

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

    % Regression matrix formulation (AR basis part)
    phiB = ones(B,K);
    for i = 1:B
        for k = 1:K
            for m = 1:M
                phiB(i,k) = phiB(i,k)*basis{m}(indx0{1}(i,m)+1,k);
            end
        end
    end
    phiY = kron(ones(1,Nr),phiB');
    BregY = kron(ones(Tdata,1),phiY);
    
    PHI_PHI = zeros(Nr*B);
    PHI_Y = zeros(Nr*B,1);
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

    RegMat{1} = PHI_PHI;
    RegMat{2} = PHI_Y;
    RegMat{3} = reshape(Y(n+1:T,:)',K*(T-n),1);

% QR or SVD method
% !!! recommended
else
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
                Breg(i,k) = Breg(i,k)*basis{m}(indx0{1}(i,m)+1,k);
            end
        end
    end
    BregY = kron(ones(1,Nr),Breg');
    BregY = kron(ones(T-n,1),BregY);

    RegMat{1} = Yreg.*BregY;                        % PHI
    RegMat{2} = reshape(Y(n+1:T,:)',K*(T-n),1);     % Ystack
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             Fitness function ('basis' search)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ftnfnc = gaPCNARXff(chromo,RegMat)

%   gaPCNARXff Computes BIC criterion for sparse Polynomial Chaos basis of a
%   PC-NARX model.
%
%   [FTNFNC] = gaPCNARXff(CHROMO,PHI)
%
%   CHROMO  : Chromosome
%   PHI     : Regression matrix of the maximum possible dimension
%
%   FTNFNC : Fitness Function to be minimized

%   Copyright 1980-2012, The MinaS Inc.
%   $ Version: 1.0 $ $Date: 26/11/2012 $

global T K Nr n B maxsize

% Initializations
PopSize = size(chromo,1);
ftnfnc = zeros(PopSize,1);

if K*(T-n)*(Nr*B) > maxsize
    PHI_PHI = RegMat{1};
    PHI_Y = RegMat{2};
    Ystack = RegMat{3};
    SSS = norm(Ystack)^2;
    for ii = 1:PopSize
        GAindex = find(kron(ones(1,Nr),chromo(ii,1:B)));
        % Parameters matrix computation
        THETA = ((PHI_PHI(GAindex,GAindex)./(T-n))\(PHI_Y(GAindex,:)./(T-n)));
        rss = SSS - 2*(PHI_Y(GAindex,:).'*THETA) + THETA'*PHI_PHI(GAindex,GAindex)*THETA;
        ftnfnc(ii,1) = log(rss/((T-n)*K)) + length(GAindex)*log(K*(T-n))/(K*(T-n));
    end
else
    % PHI = RegMat{1};
    % Ystack = RegMat{2};
    for ii = 1:PopSize
        GAindex = find(kron(ones(1,Nr),chromo(ii,1:B)));
        % AR parameter vector estimation
        THETA = RegMat{1}(:,GAindex)\RegMat{2};
        % Residuals computation
        E = RegMat{2} - RegMat{1}(:,GAindex)*THETA;
        % ftnfnc(ii,1) = 100*norm(E)^2/SSS;
        ftnfnc(ii,1) = log(var(E)) + length(GAindex)*log(K*(T-n))/(K*(T-n));
    end
end