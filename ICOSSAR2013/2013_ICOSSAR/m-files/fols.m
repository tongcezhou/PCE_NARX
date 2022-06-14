function [x, ind, r_act, A_o, x_o] = fols(A,b,r,threshold)
% FOLS Forward Orthogonal Least Quares.
%
% [x, ind] = FOLS(A,b,r) gives the solution to the least squares problem
% using only the best r regressors chosen from the ones present in matrix
% A. This function also returns in the vector ind the indexes of the
% best r regressors (i.e., the best columns of A to use).
%
% If r is equal to n, the solution x given by OLS is the same as the solution
% given by A\b, but in ind we still have the regressors sorted by their
% importance. This means that one can perform a feature selection by taking
% the first k entries in ind (k<r).
%
% THEORETICAL BACKGROUND
% Let us consider the Linear Least Squares Problem:
%
%      min   (A*x-b)'*(A*x-b)
%
% where A is a m-by-n matrix, b an m-dimentional vector and x, the unknown,
% an n-dimentional vector. The problem is well posed when m>n.
% A can be viewed as a set of columns, usually called regressors, while b
% typically is the vector of desired outputs. The solution to this problem
% can be computed in Matlab through the function MLDIVIDE as follows:
%
%   x = A\b
%
% Although this function is very powerful, it does not give any information
% on which regressor is better than others. On the other hand, Orthogonal
% Least Squares (OLS) can extract this information. It is mainly based on
% the Gram-Schmidt orthogonalization algorithm.
%
% In pattern recognition problems, usually a training matrix A_tr is given,
% associated with the vector of desired outputs b_tr, plus a test set pair
% (A_ts, b_ts) where the system has to be tested. In such cases the OLS
% function can be used as follows, with a value of r chosen by the user:
%
% [x, ind] = OLS(A_tr, b_tr, r);
%
% err_tr = (A_tr*x-b_tr)' *(A_tr*x-b_tr)/length(b_tr); % Mean Squared Error on training set.
% err_ts = (A_ts*x-b_ts)' *(A_ts*x-b_ts)/length(b_ts); % Mean Squared Error on test set.
% While a large value for r will probably give a low error on the training set,
% it could happen that a lower value of it would give a lower error on the test set.
% In such cases, the value of r which minimizes the error on the test set should be chosen.
%
% ALTERNATIVE SYNTAX
%
% [x, ind, r_act] = OLS(A,b,[], threshold) is another way to call the OLS function.
% This call requires a threshold in place of the number of regressors r.
% The threshold, which has to be a decimal number in the interval [0,1],
% indicates the maximum relative error allowed between the approximated output b_r
% and the desired output b, when the number of regressors used is r.
% The provided output r_act gives the number of regressors actually used, which
% are necessary to reduce the error untill the desired value of the threshold.
% Please not that r_act is also equal to the length of the vector ind.
% The pseudocode for this call is the following:
%
% for r=1:n
%     x_r = OLS(A,b,r);
%     b_r = A*x_r;
%     ERR = 1-sum(b_r.^2)/sum(b.^2);
%     if threshold < ERR,
% 	      break;
%     end
% end
% r_act = r;
%
% When threshold is near to 1, r tends to be 1 (a regressor tends to be sufficient).
% When threshold is near to 0, r tends to be n (all regressor tends to be needed to reach the desired accuracy). The parameter threshold is an indirect way to choose the parameter r.
%
% ALTERNATIVE SYNTAX
%
% OLS(A,b); using this syntax the OLS(A,b,r) function is called for r=1:n and the Error Reduction Ratio
% is collected and then plotted. It gives an idea on the possible values
% for threshold.
%
% ALTERNATIVE SYNTAX
%
% The function OLS computes new, orthogonal regressors in a transformed space. The new regressors
% can be extracted by using the following syntax:
%
% [x, ind, r_act, A_o] = OLS(A,b,[], threshold);
% or the following one:
%
% [x, ind, r_act, A_o] = OLS(A,b,r);
%
% The output A_o is an m-by-r_act matrix.
%
% To test the orthogonality of new generated regressors use the following instructions:
% A_o(:,u)'* A_o(:,v)
% u,v being integer values in {1,r_act}
%
% ALTERNATIVE SYNTAX
%
% [x, ind, r_act, A_o, x_o] = OLS(A,b,r); % This syntax can be used to extract the vector x_o,
% which is the vector of solutions in the orthogonally transformed space.
% This means that the length of vector x_o is r_act, as the number of columns of A_o.
% Numerically, A*x is equal to A_o*x_o.
%
% FINAL REMARK
%
% The implementation of this function follows the notation given in [1].
%
% REFERENCES
%
% [1] L. Wang, R. Langari, "Building Sugeno-yype models using fuzzy
%     discretization and Orthogonal Parameter Estimation Techniques",
%     IEEE Trans. Fuzzy Systems, vol. 3, no. 4, pp. 454-458, 1995
% [2] S.A. Billings, M. Koremberg, S. Chen, "Identification of non-linear
%     output-affine systems using an orthogonal least-squares algorithm,"
%     Int. J. Syst. Sci., vol. 19, pp. 1559-1568, 1988
% [3] M. Koremberg, S.A. Billings, Y.P. Liu, P.J. McIlroy, "Orthogonal
%     parameter estimation algorithm for nonlinear stochastic systems,"
%     Int. J. of Control, vol. 48, pp. 193-210, 1988
%
%   See also MLDIVIDE, LSQLIN, PINV, SVD, QR.

%   OLS, $Version: 0.81, August 2007
%   Author: Marco Cococcioni
%   Please report any bug to m.cococcioni <at> gmail.com

% Regression matrix size (m: number of data,  n: number of parameters)
[m,n] = size(A);

if m <= n
    error ('OLS error: The number of observations shall be strictly greater than the number of variables');
end

if nargin == 2
    SSE = zeros(n,1);
    for r = 1:n
        [x{r}, ind{r}] = fols(A,b,r); % ind{r} contains the index of the r most important regressors
        SSE(r) = (A*x{r}-b)'*(A*x{r}-b);
        b_r{r} = A*x{r};
        ERR(r) = 1-sum(b_r{r}.^2)/sum(b.^2);
    end
    
    plot(ERR,'.-b');
    ylabel('Error Reduction Ratio (ERR)');
    xlabel('Index of the most important regression variables');
    title(sprintf('Matrix A has %d rows (observations) and %d columns (regression variables)', m,n));
    set(gca,'xtick',1:n);
    set(gca,'xticklabel',ind{end});
    
    x = x{end};
    ind = ind{end};
    return
end

A = A';
cum_error = 0;
ind = [];

error_break = false;

if isempty(r),
    r = n; % In case the number of orthogonal regressors is not provided,
    % it will be set equal to the number of variables (regressors)
    error_break = true;
    % this case requires threshold
    if nargin < 4 || isempty(threshold)
        error('The threshold parameter shall be provided when r is empty');
    end
else
    if nargin > 3 && ~isempty(threshold)
        error('If r is provided, the threshold parameter has to be empty');
    end
end

if isa(r,'integer') || (r > n) || (r <= 0)
    error('The number of regressors r must be a positive integer in the range [1,size(A,2)]');
end

% Series Sum of Squares
SSS = norm(b)^2;

A_o = zeros(n,m);
alpha_tot = zeros(n,n);
x_o = zeros(1,n);

% Gram-Schmidt orthogonalization -- Eq. (24)
% Selection of the first term
wi = A;
wi2 = sum(wi.^2,2);
gi = wi*b./wi2;
erri = (gi.^2).*wi2/SSS;
% The term producing the maximum ERRi is selected
[value, index] = max(erri);
cum_error = cum_error + value;
reg_err = 1 - cum_error;
ind = [ind index];
avail_index = [1:(index-1) (index+1):n];
A_o(1,:) = A(index,:);
x_o(1,1) = gi(index);

if (error_break && (reg_err <= threshold))
    r = 1;   % avoid entering next cycle
end

r_act = 1;

for m = 2:r
    kk = 0;
    sav = size(avail_index,2);
    alpha = zeros(n,sav);
    
    for av = avail_index
        kk = kk+1;
        for i = 1:m-1
            alpha(i,kk) = A_o(i,:)*(A(av,:).')/sum(A_o(i,:).^2);
        end
    end
    
    wi = A(avail_index,:)-(alpha(1:m-1,:)')*A_o(1:m-1,:);
    wi2 = sum(wi.^2,2);
    gi = wi*b./wi2;
    erri = (gi.^2).*wi2/SSS;
    [value,index] = max(erri);
    
    cum_error = cum_error + value;
    reg_err = 1 - cum_error;
    ind(m) = avail_index(1,index);
    davail_index = size(avail_index,2);
    avail_index = [avail_index(1:(index-1)) avail_index((index+1):davail_index)];
    A_o(m,:) = wi(index,:);
    x_o(1,m) = gi(index);
    alpha_tot(:,m) = alpha(:,index);
    r_act = m;
    if (error_break && (reg_err <= threshold))
        break
    end
end

sav = size(ind,2);
x = zeros(1,n);
x(1,sav) = x_o(1,sav);                                          % Eq. (21)
for i = sav-1:-1:1
    x(1,i) = x_o(1,i)-alpha_tot(i,i+1:sav)*(x(1,i+1:sav).');    % Eq. (22)
end

tmp = x;
x = zeros(1,n);
for i = 1:sav
    x(1,ind(i)) = tmp(1,i);
end

x = x';
A_o = A_o';
A_o = A_o(:,1:r_act);
x_o = x_o(1,1:r_act)';
end
