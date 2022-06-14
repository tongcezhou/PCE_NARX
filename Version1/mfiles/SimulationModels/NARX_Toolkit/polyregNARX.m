function regstr = polyregNARX(orders,degree,q,maxterms,addterms)

%   polyregNARX forms the polynomial nonlinear regressors f(y[t-j], x[t-k])
%       for an NARX model:
%
%       y[t] + \sum_{i=1}^{Nr} theta_i * f(y[t-j], x[t-k]) = e[t]
%
%       The regressors constructed are polynomials of maximum total degree
%       equal to P, maximum quasi-norm Q, and maximum number of terms in
%       each regressor equal to MAXTERMS.
%
%	REGSTR = polyregNARX(ORDERS,P,Q,MAXTERMS)
%
%   ORDERS = [Na Nd Nb] : The AR/X orders of the NARX input-output model
%   P : Maximum total polynomial degree for each regressor, e.g. for
%               '(Y(t-1,:).^i)*(X(t,:).^j)'  --->  i+j < P
%   Q : quasi-norms order (help combinations)
%   MAXTERMS : Maximum number of terms included in each regressor
%
%   REGSTR : Cell array consisting of the constructed regressors

%   Copyright 1980-2013, The MinaS Inc.
%   $ Version: 1.0 $ $Date: 18/04/2013 $

% Check input variables
if isempty(maxterms),      maxterms = 1;   end
if isempty(q),             q = 1;          end


% % Check if additional terms should be incorporated
% if nargin == 5
%     maxterms = 1;
% end

% NAR case
% =================
if length(orders) == 1
    % AR order
    na = orders(1);
    if maxterms == 1
        POWindx = [zeros(1,na) ; kron([1:degree]',eye(na))];
    else
        POWindx = combinations(na,degree,q);
        POWindx = POWindx(sum(logical(POWindx),2) <= maxterms,:);
    end
    
    regstr = cell(size(POWindx,1),1);
    regstr{1} = 'ones(length(t),size(Y,2))';
    for j = 2:size(POWindx,1)
        for k = 1:na
            if POWindx(j,k)>0
                if isempty(regstr{j})
                    regstr{j} = [regstr{j},'(Y(t-',num2str(k),',:).^',num2str(POWindx(j,k)),')'];
                else
                    regstr{j} = [regstr{j},'.*(Y(t-',num2str(k),',:).^',num2str(POWindx(j,k)),')'];
                end
            end
        end
    end
    % NARX case
    % =================
else
    % AR order
    na = orders(1);
    % Input delay order
    nd = orders(2);
    % X order
    nb = orders(3);
    
    if maxterms == 1
        POWindx = [zeros(1,(na+nb-nd+1)) ; kron([1:degree]',eye(na+nb-nd+1))];
    else
        POWindx = combinations(na+nb-nd+1,degree,q);
        POWindx = POWindx(sum(logical(POWindx),2) <= maxterms,:);
    end
    regstr = cell(size(POWindx,1),1);
    regstr{1} = 'ones(length(t),size(Y,2))';
    for j = 2:size(POWindx,1)
        for k = 1:na+nb-nd+1
            if POWindx(j,k)>0
                if k <= na
                    if isempty(regstr{j})
                        regstr{j} = [regstr{j},'(Y(t-',num2str(k),',:).^',num2str(POWindx(j,k)),')'];
                    else
                        regstr{j} = [regstr{j},'.*(Y(t-',num2str(k),',:).^',num2str(POWindx(j,k)),')'];
                    end
                else
                    if isempty(regstr{j})
                        regstr{j} = [regstr{j},'(X(t-',num2str(k-na+nd-1),',:).^',num2str(POWindx(j,k)),')'];
                    else
                        regstr{j} = [regstr{j},'.*(X(t-',num2str(k-na+nd-1),',:).^',num2str(POWindx(j,k)),')'];
                    end
                end
            end
        end
    end
end



% Check if additional terms should be incorporated
if nargin == 5
    curlen = size(regstr,1);
    if iscell(addterms)
        p = size(regstr,1);
        for i = 2:curlen
            for j = 1:length(addterms)
                p = p+1;
                regstr{p} = [regstr{i},'.*',addterms{j}];
            end
        end
    elseif strcmp(addterms,'abs')
        p = size(regstr,1);
        for i = 2:curlen
            for j = 1:na+nb-nd+1
                p = p+1;
                if k <= na
                    regstr{p} = [regstr{i},'.*','abs(Y(t-',num2str(k),',:))'];
                else
                    regstr{p} = [regstr{i},'.*','abs(X(t-',num2str(k-na+nd-1),',:))'];
                end
            end
        end
    end
end