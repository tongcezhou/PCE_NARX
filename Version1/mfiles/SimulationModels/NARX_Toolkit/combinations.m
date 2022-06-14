function [indx,P,J] = combinations(M,p,q) 

%   COMBINATIONS Computes the vector integer index containing the degree of
%   polynomial basis functions included in each term.
%
%   [INDX,P,J] = COMBINATIONS(M,P,Q)
%
%   M : Number of variables
%   P : Total maximum polynomial degree
%   Q : quasi-norms order [1]
%             
%   INDX : Index of the polynomial basis
%   P : Index total polynomial degree
%   J : Index interaction order 
%
%   [1] Blatman, G. & Sudret, B. Adaptive Sparse Polynomial Chaos Expansion 
%       Based on Least Angle Regression Journal of Computational Physics,
%       2011, Vol. 230, pp. 2345-2367.

%   Copyright 1980-2013, The MinaS Inc.
%   $ Version: 1.01 $ $Date: 18/04/2013 $

if nargin<3
    q = 1;
end

% Current row of the indx matrix
row_cur = zeros(M,1);
for ii=1:M
    row_cur(ii)= round(factorial(ii+p)/(factorial(ii)*factorial(p)));
end

indx = zeros(row_cur(end),M);
indx(1:p+1,1) = (0:p)';

for ii = 2:M
    indx_cur = row_cur(ii-1)+1;
    for jj = 0:p-1
        J = find(sum(indx(1:row_cur(ii-1),:),2) == jj);
        for m = 1:p-jj
            indx(indx_cur:indx_cur + length(J)-1,:) = indx(J,:);
            indx(indx_cur:indx_cur + length(J)-1,ii) = m ;
            indx_cur = indx_cur + length(J);
        end
    end
end

if q < 1
    indx_in = [];
    for ii = 1:size(indx,1)
        if norm(indx(ii,:),q) <= p
            indx_in = [indx_in ; ii];
        end
    end
    indx = indx(indx_in,:);
end

[~,jj] = sort(sum(indx,2));
indx = indx(jj,:);

[~,jj] = sort(sum(indx>0,2));
indx = indx(jj,:);

P = sum(indx,2);
J = sum(indx>0,2);