function [SEstructure,IndxRem,IndxRmv] = SEreduction(Y,X,orders,regressors,options)

%   BICreduction calculates the BIC criterion for the second stage of the
%   heuristic scheme for FS-TARMA model structure selection. The Functional
%   Series AutoRegressive Moving Average model is:
%
%   y[t]+\sum_{i=1}^{na} a_i[t]*y[t-i]=e[t]+\sum_{i=1}^{nc} c_i[t]*e[t-i]
%
%   [BICstruct,INDnew] = BICreduction(Y,ORDERS,IND,options,InitRes)
%
%   Y : Output data (1 x N)
%   ORDERS : The FS-TARMA output model orders [Na Nc] (1 x 2)
%           Na : The AR order
%           Nc : The MA order
%   IND : basis function indices (structure)
%           Ba : AR basis function indices vector
%           Bc : MA basis function indices vector
%           Bs : Variance basis function indices vector
%   OPTIONS : estimation method options (structure)
%             basis : 'cheby' or 'sinus' or '-haar'
%             est_meth : 'sum' or '-qr' or 'svd' or '-ml', (prefer -qr if possible)
%             lossfun : 'ols' (ordinary least squares) or '-ml' (maximum likelihood) or 'wls' (weighted least squares)
%             s_est : 'np' (non-parametric)  or '-p' (parametric) or '-g' (std estimation)
%             s_ref : 'y' or 'n' (s ML refinement)
%             res_t0 : number of initial resiuals excluded from the criteria computation
%             t0 : initial time instant
%             maxsize : maximum number of elements allowed in an array ('sum' estimation method)
%             M : window length = 2*M+1
%             s_ml_options : optimization OPTIONS structure
%             ml_options : optimization OPTIONS structure
%             parameters_ref : parameter vector to be optimized ('-ar' or '-ma' or 'all')
%   INITRES : Initial residuals (2SLS and PA starting from 2nd stage)
%
%   BICstruct : structure
%       VALUES : BIC values obtained during basis removal process;
%       Basis : sequence of strings indicating from which subspace a
%               basis was removed ('a'/'c'/'s' for AR/MA/variance)
%       Removed : Sequence with the indeces of the basis removed
%   INDnew : finally selected basis function indices (structure)
%           Ba : AR basis function indices vector
%           Bc : MA basis function indices vector
%           Bs : Variance basis function indices vector
%
%   Copyright 1980-2011, The MinaS Inc.
%   $ Version: 1.1 $ $Date: 15/11/2011 $

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isfield(options,'tolfun') && ~isfield(options,'Nr')
    options.Nr = 10;
end
if ~isfield(options,'focus'),   options.focus = 'prediction';   end
%======================================================================
% NARX CASE
%======================================================================
% AR order
na = orders(1);
% Input delay order
nd = orders(2);
% X order
nb = orders(3);
% Initial time instant
n = max(na,nb);
% Number of samples
T = size(Y,1);

% Initial regressor vector
REG0 = regressors;
% Sum os squared series
SSS = norm(Y)^2;
% Initial FS-TAR model estimation
th = narx(Y,X,[na nd nb],regressors,[],[],[],[],options);
Ysim = narxsim(X,Y(1:na,:),na,regressors,th);

% RSS/SSS criterion
if strcmp(options.criterion,'rss')
    SE0 = 100*norm(Y-Ysim)^2/SSS;
% MNSE
else
    SE0 = mean(abs(Y-Ysim)./(1+abs(Y)));    
end

if isnan(SE0)
    SE0 = inf;
end
% minSE = criteria.rss_sss;
minSE = SE0;
% Iterations until no further reduction of the BIC value is possible
% (AR functional subspace reduction)
SEindx = [];
SEstructure = [];
SEstructure.initvalue = SE0;

Ymeas = Y;
Xmeas = X;

if isfield(options,'tolfun')
    k = 1; % counter
    while minSE<SE0 || abs(minSE-SE0)/(1+SE0) < options.tolfun
        if k > 1
            regressors = regressors([1:SEindx-1 SEindx+1:end]);
        end
        % Number of regressors
        Nr = size(regressors,1);
        clear Yreg Y RegMat X
        Y = Ymeas;
        X = Xmeas;
        t = n+1:T;
        % Y regressor part
        Yreg = zeros(T-n,Nr);
        for j = 1:Nr
            Yreg(1:T-n,j) = eval(regressors{j});
        end
        RegMat{1} = Yreg;         	% PHI
        RegMat{2} = Y(n+1:T,:);    	% Ystack
        
        regrowstr = [];
        % Concentrate regressors in a single string
        regrowstr = strcat('[',regressors{1});
        for j = 2:size(regressors,1)
            regrowstr = strcat(regrowstr,' ; ',regressors{j});
        end
        regrowstr = strcat(regrowstr,']');
        
        pa = length(regressors);
        THETAindexed = zeros(pa,pa);
        % Initial model
        SE = inf*ones(pa,1);   % Setting every element equal to inf
        for ii = 1:pa
            aux= ones(pa,1);
            aux(ii) = 0;
            indx = find(aux == 1);
            % AR parameter vector estimation
            THETA = RegMat{1}(:,indx)\RegMat{2};
            % Indexed complete theta vector
            THETAindexed(indx,ii) = THETA;
        end
        
        % Create copies of X
        X =  kron(ones(1,pa),Xmeas);
        % Initialize set of simulated responses
        Y = zeros(T,pa);
        Y(1:n,:) = kron(ones(1,pa),Ymeas(1:n,1));
        % Simulation
        for t = n+1:T
            Y(t,:) = sum(eval(regrowstr).*THETAindexed);
        end
        
        % Criteria computation
        for ii = 1:pa
            E = Ymeas - Y(:,ii);
            if isnan(E(end)) || isinf(E(end))
                SE(ii,1) = nan;
            else
                % RSS/SSS criterion
                if strcmp(options.criterion,'rss')
                    SE(ii,1) = 100*norm(E)^2/SSS;
                % MNSE
                else
                    SE(ii,1) = mean(abs(E)./(1+abs(Ymeas)));
                end
             end
        end
        [minSE,SEindx] = min(SE);
        % If minBIC is smaller than the previous value then the
        % corresponding basis will be removed
        if ~isnan(minSE) && ~isinf(minSE) && (minSE<SE0 || abs(minSE-SE0)/(1+SE0) < options.tolfun)
            SEstructure.values(k) = minSE;
            SEstructure.removed(k) = regressors(SEindx);
            k = k+1;
        else
            break;
        end
    end
else
    IndxRmv = zeros(options.Nr,1);
    k = 1; % counter
    while k <= options.Nr
        if k > 1
            regressors = regressors([1:SEindx-1 SEindx+1:end]);
        end
        disp(k);
        % Number of regressors
        Nr = size(regressors,1);
        clear Yreg Y RegMat X
        Y = Ymeas;
        X = Xmeas;
        t = n+1:T;
        % Y regressor part
        Yreg = zeros(T-n,Nr);
        for j = 1:Nr
            Yreg(1:T-n,j) = eval(regressors{j});
        end
        RegMat{1} = Yreg;         	% PHI
        RegMat{2} = Y(n+1:T,:);    	% Ystack
        
        regrowstr = [];
        % Concentrate regressors in a single string
        regrowstr = strcat('[',regressors{1});
        for j = 2:size(regressors,1)
            regrowstr = strcat(regrowstr,' ; ',regressors{j});
        end
        regrowstr = strcat(regrowstr,']');
        
        pa = length(regressors);
        THETAindexed = zeros(pa,pa);
        % Initial model
        SE = inf*ones(pa,1);   % Setting every element equal to inf
        for ii = 1:pa
            aux = ones(pa,1);
            aux(ii) = 0;
            indx = find(aux == 1);
            % AR parameter vector estimation
            THETA = RegMat{1}(:,indx)\RegMat{2};
            % Indexed complete theta vector
            THETAindexed(indx,ii) = THETA;
        end
        
        % Create copies of X
        X =  kron(ones(1,pa),Xmeas);
        % Initialize set of simulated responses
        Y = zeros(T,pa);
        Y(1:n,:) = kron(ones(1,pa),Ymeas(1:n,1));
        % Simulation
        for t = n+1:T
            Y(t,:) = sum(eval(regrowstr).*THETAindexed);
        end
        
        % Criteria computation
        for ii = 1:pa
            E = Ymeas - Y(:,ii);
            if isnan(E(end)) || isinf(E(end))
                SE(ii,1) = nan;
            else
                % RSS/SSS criterion
                if strcmp(options.criterion,'rss')
                    SE(ii,1) = 100*norm(E)^2/SSS;
                % MNSE
                else
                    SE(ii,1) = mean(abs(E)./(1+abs(Ymeas)));
                end
            end
        end
        [minSE,SEindx] = min(SE);
        % If minBIC is smaller than the previous value then the
        % corresponding basis will be removed
        if ~isnan(minSE) && ~isinf(minSE) && (k <= options.Nr)
            SEstructure.values(k) = minSE;
            SEstructure.removed(k) = regressors(SEindx);
            for q = 1:size(REG0,1)
                if strcmp(REG0(q),regressors(SEindx))
                    IndxRmv(k) = q;
                    break;
                end
            end
            k = k+1;
        else
            break;
        end
    end
end

IndxRem = zeros(size(REG0,1),1);
for j = 1:size(REG0,1)
    for k = 1:size(regressors,1)
        if strcmp(REG0(j),regressors(k))
            IndxRem(j) = 1;
            break;
        end
    end
end

% %======================================================================
% % NARX CASE
% %======================================================================
% % AR order
% na = orders(1);
% % Input delay order
% nd = orders(2);
% % X order
% nb = orders(3);
% % Initial time instant
% n = max(na,nb);
%
% SSS = norm(Y)^2;
%
% % Initial FS-TAR model estimation
% SE0 = inf;
% % minSE = criteria.rss_sss;
% minSE = 100;
% % Iterations until no further reduction of the BIC value is possible
% % (AR functional subspace reduction)
% CURreg = [];
% SEindx = [];
%
% k = 1; % counter
% while ~isempty(regressors)%(minSE < SE0)
%     if k > 1
%         SE0 = SE(SEindx);
%         CURreg = [CURreg ; regressors(SEindx)];
%         regressors = regressors([1:SEindx-1 SEindx+1:end]);
%     end
%     pa = length(regressors);
%     % Initial model
%     SE = inf*ones(pa,1);   % Setting every element equal to inf
%     for ii = 1:pa
%         % Remove from the AR functional subspace one basis at a time
%         reg_cur = [CURreg ; regressors(ii)];
%         % Estimated the new NARX model
%         [th,~,criteria] = narx(Y,X,[na nd nb],reg_cur,[],[],[],[],options);
%         Ysim = narxsim(X,Y(1:na,:),na,reg_cur,th);
%         if ~isnan(norm(Ysim)) && ~isinf(norm(Ysim))
%             % Calculating the BIC parameters
%             %SE(ii) = criteria.rss_sss;
%             SE(ii) = 100*norm(Ysim)^2/SSS;
%         end
%         % Calculating the APD criterion
%     end
%     [minSE,SEindx] = min(SE);
%     % If minBIC is smaller than the previous value then the
%     % corresponding basis will be removed
%     if minSE < 100
%         SEstructure.values(k) = minSE;
%         SEstructure.added(k) = regressors(SEindx);
%         disp(minSE)
%         k = k+1;
%     end
% end