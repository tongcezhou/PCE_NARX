function regstr = polyregNARX(orders,degree,q,maxterms)

%   polyregNARX forms the polynomial nonlinear regressors f(y[t-j], x[t-k])
%       for an NARX model:
%
%       y[t] + \sum_{i=1}^{Nr} theta_i * f(y[t-j], x[t-k]) = e[t]
%
%       The regressors constructed are polynomials of maximum total degree
%       equal to P, maximum quasi-norm Q, and maximum number of terms in
%       each regressor equal to MAXTERMS.
%
%	REGSTR = polyregNARX(Y,X,ORDERS,P,Q,MAXTERMS)
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

% NAR case
% =================
if length(orders) == 1
    % AR order
    na = orders(1);
    POWindx = combinations(na,degree,q);
    if ~isempty(maxterms)
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
    
    POWindx = combinations(na+nb-nd+1,degree,q);
    if ~isempty(maxterms)
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