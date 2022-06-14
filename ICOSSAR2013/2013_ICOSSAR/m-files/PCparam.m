function [phi,an] = PCparam(Q,p,indx,A,options)

%   PCEparam Computes the parameter d_i(x) of an estimated Polynomial 
%   Chaos model
%
%   [PHI,An] = PCEparam(Q,P,IND,A,OPTIONS)
%
%   Q : Input data (M x K)
%   P : Maximum basis degree
%   IND : basis function indices (B x M)
%   A : Estrimated PC basis coefficients of projection vector (n*B x 1)
%   OPTIONS : (structure)
%             basis : 'hermi'/'cheby'/'legen'...
%
%   Phi : Basis matrix (B x K)
%   An : PC coefficients (K x n*B)
              
%   Copyright 1980-2013, The MinaS Inc.
%   $ Version: 1.0 $ $Date: 19/04/2013 $

% Number of variables
M = size(Q,1);
% Number of experiments
K = size(Q,2);

% Basis constraction: size(basis_i) = p x K
basis = cell(M,1);
for i = 1:M
    basis{i} = PCbasis(Q(i,:),p,options);
end

B = size(indx,1);
% Regression matrix formulation (basis part)
phi = ones(B,K);
for i = 1:B
    for k = 1:K
        for m = 1:M
            phi(i,k) = phi(i,k)*basis{m}(indx(i,m)+1,k);
        end
    end
end

% Number of parameters 
n = length(A)/B;
% Parameters computation
if ~isempty(A) && nargout > 1 
    an = phi'*reshape(A,B,n);
end