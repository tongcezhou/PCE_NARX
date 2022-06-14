function [GAreg,indx,RES,criteria,GAoutput] = gaNARX(Y,X,orders,regressors,options)

% [GAreg,INDX,E,CRITERIA,GAoutput] = gaNARX(Y,X,ORDERS,REGRESSORS,OPTIONS)
%
%   gaPCARX: GA-based model structures selection for a Nonlinear
%   AutoRegressive with eXogenous input model:
%
%   y[t] + \sum_{i=1}^{Nr} theta_i * f(y[t-j], x[t-k]) = e[t]
%
%   Y : Output data (T X 1)
%   X : Excitation input data (T X 1)
%   ORDERS = [Na Nd Nb] : The AR/X orders of the PC-ARX input-output model
%   REGRESSORS : Cell array consisting of the regressors in string format e.g.
%       '(Y(t-2,:).^3)*(X(t,:))'
%   OPTIONS : estimation method options (structure)
%       focus : 'prediction'/'simulation'
%       criterion : 'bic'/'rss'/'mnse'
%   OPTIONS.GA : Genetic algorithm options (structure)
%
%   GAreg : Cell array consisting of the finally selected regressors
%   INDX : The finally selected regressors index
%   E : estimated residual series (T x 1)
%   CRITERIA :
%    	SPP : Samples Per Parameter
%     	CONNUM : Regressor matrix condition number
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
%
%   * If REGSTR is empty the regressor terms are formed by polynomial basis
%          functions of maximum total degree equal to 3, maximum quasi-norm 1,
%          and maximum number of terms in each regressor equal to 1.

%   Copyright 1980-2014, The MinaS Inc.
%   $ Version: 1.25 $ $Date: 07/01/2014 $

% Global variables
global T n na nb nd

if strcmp(options.warnings,'off')
    % Turning off warning for rank deficiency
    warning('off','MATLAB:rankDeficientMatrix');
    warning('off','MATLAB:nearlySingularMatrix');
end

if ~(isfield(options,'focus')),                  options.focus = 'prediction' ;             end
% Genetic Algorithm options
if ~(isfield(options.GA,'PopulationSize')),      options.GA.PopulationSize = 100 ;          end
if ~(isfield(options.GA,'Generations')),         options.GA.Generations = 200;              end
if ~(isfield(options.GA,'CrossoverFraction')),   options.GA.CrossoverFraction = 0.7;        end
if ~(isfield(options.GA,'Display')),             options.GA.Display = 'iter';               end
if ~(isfield(options.GA,'PlotFcns')),            options.GA.PlotFcns = {@gaplotbestf,@gaplotbestindiv};  end
if ~(isfield(options.GA,'StallGenLimit')),       options.GA.StallGenLimit = 20;             end
if ~(isfield(options.GA,'EliteCount')),
    options.GA.EliteCount = round(0.05*options.GA.PopulationSize);                          end
if ~(isfield(options.GA,'MutationFcn')),
    options.GA.MutationFcn = {@mutationuniform,0.3};                                        end

options.GA.PopulationType = 'bitstring';
options.GA.Vectorized = 'on';
options.GA.StallTimeLimit = Inf;
GAoptions = gaoptimset(options.GA);

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
if isempty(regressors)
    regressors = polyregNARX([na nd nb],3,1,1);
end
% Initial regressors
customreg = regressors;
% Number of regressors
Nr = length(regressors);


t = n+1:T;
% Y regressor part
Yreg = zeros(T-n,Nr);
for j = 1:Nr
    Yreg(1:T-n,j) = eval(customreg{j});
end
RegMat{1} = Yreg;         	% PHI
RegMat{2} = Y(n+1:T,:);    	% Ystack

% Gene length
gene = size(RegMat{1},2);
% Fitness function computation
% PARFOR loop
if strcmp(options.parfor,'y')
    matlabpool('open',4)
end
[indiv,GAoutput.fval,GAoutput.reason,GAoutput.output,GAoutput.population,GAoutput.scores] = ...
    ga(@(gene) gaNARXff(gene,Y,X,RegMat,customreg,options),gene,[],[],[],[],[],[],[],GAoptions);

if strcmp(options.parfor,'y')
    matlabpool('close')
end


indx = find(indiv == 1);
Nr = length(indx);
GAreg = cell(Nr,1);
for j = 1:Nr
    GAreg{j} = regressors{indx(j)};
end

THETA = RegMat{1}(:,indx)\RegMat{2};
% Simulated response computation
Ysim = narxsim(X,Y(1:na),na,GAreg,THETA);
% Residual computation
RES = Y-Ysim;

% Criteria Computation
%%%%%%%%%%%%%%%%%%%%%%
SSS = norm(Y(n+1:T))^2;
est_param = length(THETA);
criteria.spp = (T-n)/est_param;
criteria.rss = norm(RES(n+1:T))^2;
criteria.rss_sss = 100*criteria.rss/SSS;
criteria.mse = criteria.rss/(T-n);
criteria.mnse = mean(abs(RES(n+1:T))./(1+abs(Y(n+1:T))));
criteria.bic = log(var(RES)) + est_param*log(T-n)/(T-n);

% Turning on warning for rank deficiency
warning('on','MATLAB:rankDeficientMatrix');
warning('on','MATLAB:nearlySingularMatrix');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             Fitness function ('basis' search)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ftnfnc = gaNARXff(chromo,Ymeas,X,RegMat,regressors,options)

%   gaNARXff Computes the fitness criterion for the selection of the optimal regressors

global T n

% Initializations
PopSize = size(chromo,1);
ftnfnc = zeros(PopSize,1);

% PHI = RegMat{1};
% Ystack = RegMat{2};

% Number of samples
Nos = length(RegMat{2});

% PARFOR loop
if strcmp(options.parfor,'y')
    % Prediction-Error criterion
    if strcmp(options.focus,'prediction')
        % Series Sum of Squares
        SSS = norm(Ymeas)^2;
        
        parfor ii = 1:PopSize
            indx = find(chromo(ii,:)==1);
            % AR parameter vector estimation
            THETA = RegMat{1}(:,indx)\RegMat{2};
            % Residuals computation
            E = RegMat{2} - RegMat{1}(:,indx)*THETA;
            
            % BIC criterion
            if strcmp(options.criterion,'bic')
                ftnfnc(ii,1) = log(var(E)) + length(indx)*log(Nos)/(Nos);
                % RSS/SSS criterion
            elseif strcmp(options.criterion,'rss')
                ftnfnc(ii,1) = 100*norm(E)^2/SSS;
                % MNSE criterion
            else
                ftnfnc(ii,1) = mean(abs(E)./(1+abs(Ymeas)));
            end
        end
        % Simulation-Error criterion
    else
        % Concentrate regressors in a single string
        regrowstr = strcat('[',regressors{1});
        for j = 2:size(regressors,1)
            regrowstr = strcat(regrowstr,' ; ',regressors{j});
        end
        regrowstr = strcat(regrowstr,']');
        
        % Create theta vector full of zeros (the chromosome index will be used
        % to identify the regressors included)
        THETAindxexed = zeros(size(regressors,1),PopSize);
        
        for ii = 1:PopSize
            indx = find(chromo(ii,:)==1);
            % AR parameter vector estimation
            THETA = RegMat{1}(:,indx)\RegMat{2};
            % Indexed complete theta vector
            THETAindxexed(indx,ii) = THETA;
        end
        % Create copies of X
        X =  kron(ones(1,PopSize),X);
        % Initialize set of simulated responses
        Y = zeros(T,PopSize);
        Y(1:n,:) = kron(ones(1,PopSize),Ymeas(1:n,1));
        % Simulation
        for t = n+1:T
            Y(t,:) = sum(eval(regrowstr).*THETAindxexed);
        end
        % Series Sum of Squares
        SSS = norm(Ymeas)^2;
        % Criteria computation
        for ii = 1:PopSize
            indx = find(chromo(ii,:)==1);
            E = Ymeas - Y(:,ii);
            if isnan(E(end)) || isinf(E(end))
                ftnfnc(ii,1) = NaN;
            else
                % BIC criterion
                if strcmp(options.criterion,'bic')
                    ftnfnc(ii,1) = log(var(E)) + length(indx)*log(Nos)/(Nos);
                    % RSS/SSS criterion
                elseif strcmp(options.criterion,'rss')
                    ftnfnc(ii,1) = 100*norm(E)^2/SSS;
                    % MNSE criterion
                else
                    ftnfnc(ii,1) = mean(abs(E)./(1+abs(Ymeas)));
                end
            end
        end
    end
    
% Conventional for loop
else
    
    % Prediction-Error criterion
    if strcmp(options.focus,'prediction')
        % Series Sum of Squares
        SSS = norm(Ymeas)^2;
        
        parfor ii = 1:PopSize
            indx = find(chromo(ii,:)==1);
            % AR parameter vector estimation
            THETA = RegMat{1}(:,indx)\RegMat{2};
            % Residuals computation
            E = RegMat{2} - RegMat{1}(:,indx)*THETA;
            
            % BIC criterion
            if strcmp(options.criterion,'bic')
                ftnfnc(ii,1) = log(var(E)) + length(indx)*log(Nos)/(Nos);
                % RSS/SSS criterion
            elseif strcmp(options.criterion,'rss')
                ftnfnc(ii,1) = 100*norm(E)^2/SSS;
                % MNSE criterion
            else
                ftnfnc(ii,1) = mean(abs(E)./(1+abs(Ymeas)));
            end
        end
        % Simulation-Error criterion
    else
        % Concentrate regressors in a single string
        regrowstr = strcat('[',regressors{1});
        for j = 2:size(regressors,1)
            regrowstr = strcat(regrowstr,' ; ',regressors{j});
        end
        regrowstr = strcat(regrowstr,']');
        
        % Create theta vector full of zeros (the chromosome index will be used
        % to identify the regressors included)
        THETAindxexed = zeros(size(regressors,1),PopSize);
        
        for ii = 1:PopSize
            indx = find(chromo(ii,:)==1);
            % AR parameter vector estimation
            THETA = RegMat{1}(:,indx)\RegMat{2};
            % Indexed complete theta vector
            THETAindxexed(indx,ii) = THETA;
        end
        % Create copies of X
        X =  kron(ones(1,PopSize),X);
        % Initialize set of simulated responses
        Y = zeros(T,PopSize);
        Y(1:n,:) = kron(ones(1,PopSize),Ymeas(1:n,1));
        % Simulation
        for t = n+1:T
            Y(t,:) = sum(eval(regrowstr).*THETAindxexed);
        end
        % Series Sum of Squares
        SSS = norm(Ymeas)^2;
        % Criteria computation
        for ii = 1:PopSize
            indx = find(chromo(ii,:)==1);
            E = Ymeas - Y(:,ii);
            if isnan(E(end)) || isinf(E(end))
                ftnfnc(ii,1) = NaN;
            else
                % BIC criterion
                if strcmp(options.criterion,'bic')
                    ftnfnc(ii,1) = log(var(E)) + length(indx)*log(Nos)/(Nos);
                    % RSS/SSS criterion
                elseif strcmp(options.criterion,'rss')
                    ftnfnc(ii,1) = 100*norm(E)^2/SSS;
                    % MNSE criterion
                else
                    ftnfnc(ii,1) = mean(abs(E)./(1+abs(Ymeas)));
                end
            end
        end
    end
end
