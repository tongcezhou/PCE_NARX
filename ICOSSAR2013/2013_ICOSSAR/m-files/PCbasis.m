function basis = PCbasis(Q,P,options)

%   PCbasis computes the functional basis to be used for the Polynomial 
%   Chaos (PC) expansion.
%
%   BASIS = PCbasis(Q,P,OPTIONS)
%
%   Q : Input matrix (K x M)
%   P : Maximum polynomial degree
%   OPTIONS : basis options (structure)
%             basis : 'hermi'/'lague'/'cheby'/'legen'/'jacob'
%             ortho : 'y' or 'n' (orthogonalization)
%             param  = [a b] : coefficients of the Jacobi polynomials 
%                   or [a] : coefficient of the Laguerre polynomials
%             
%   BASIS : Functional basis (P x K)
%           
% [1] C. Soize and R. Ghanem, "Physical Systems with Random Uncertainties:
%       Chaos Representations with Arbitrary Probability Measure", SIAM 
%       Journal of Scientific Computing, Vol. 26, No. 2, pp. 395-410, 2004.

%   Copyright 1980-2012, The MinaS Inc.
%   $ Version: 1.01 $ $Date: 06/07/2012 $

if nargin<3,                        options.basis = 'hermi';    end
if ~(isfield(options,'basis')),     options.basis = 'hermi';    end
if ~(isfield(options,'ortho')),     options.ortho = 'y';        end

% Functional basis computation
switch lower(options.basis)
    case 'hermi' % Hermite polynomials (continuous)
        basis = basis_hermite(Q,P,options.ortho);
    case 'lague' % Laguerre polynomials (continuous)
        if ~(isfield(options,'param'));
            options.param = 0;
        end
        alpha = options.param(1);
        basis = basis_laguerre(alpha,Q,P,options.ortho);
    case 'cheby' % Chebyshev polynomials (continuous)
        basis = basis_chebyshev(Q,P,options.ortho);
    case 'legen' % Legendre polynomials (continuous)
        basis = basis_jacobi(0,0,Q,P,options.ortho);
    case 'jacob' % Jacobi polynomials (continuous)
        if ~(isfield(options,'param'));
            options.param = [-0.3; 0.8];
        end
        alpha = options.param(1);
        beta =  options.param(2);
        basis = basis_jacobi(alpha,beta,Q,P,options.ortho);
    otherwise
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 Hermite Polynomial Basis Construction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function basis = basis_hermite(Q,P,ortho)

% Support: Real axis
% Associated Probability measure: Gaussian

% Number of experiments
K = size(Q,2);
% Basis variable preallocation
basis = zeros(P+1,K);

basis(1,:) = 1; % H_0(Q)
basis(2,:) = Q; % H_1(Q)

for i = 2:P    
    basis(i+1,:) = Q.*basis(i,:)-(i-1)*basis(i-1,:);   % H_k(Q)
end

if strcmp(ortho,'y')
    for i = 0:P
        basis(i+1,:) = basis(i+1,:)/sqrt(factorial(i));
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Laguerre Polynomial Basis Construction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function basis = basis_laguerre(a,Q,P,ortho)

% Support: ]0,+infty[
% Associated Probability measure: Uniform
% a > -1

% Number of experiments
K = size(Q,2);
% Basis variable preallocation
basis = zeros(P+1,K);

basis(1,:) = 1;     % L_0(Q)
basis(2,:) = a+1-Q; % L_1(Q)

for i = 2:P    
    basis(i+1,:) = ((2*i-1+a-Q).*basis(i,:)-(i-1+a)*basis(i-1,:))/i;   % L_k(Q)
end

if strcmp(ortho,'y')
    for i = 0:P
        basis(i+1,:) = basis(i+1,:)*sqrt(factorial(i)/gamma(a+i+1));
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Chebyshev Polynomial Basis Construction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function basis = basis_chebyshev(Q,P,ortho)

% Support: ]-1,1[
% Associated Probability measure: Chebyshev

% Number of experiments
K = size(Q,2);
% Basis variable preallocation
basis = zeros(P+1,K);

basis(1,:) = 1; % T_0(Q)
basis(2,:) = Q; % T_1(Q)

for i = 2:P    
    basis(i+1,:) = 2*Q.*basis(i,:) - basis(i-1,:);   % T_k(Q)
end

if strcmp(ortho,'y')
    for i = 0:P
        if i == 0,  delta = 1;  else    delta = 0;  end
        basis(i+1,:) = basis(i+1,:)*sqrt(2/(pi*(1 + delta)));
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Jacobi Polynomial Basis Construction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function basis = basis_jacobi(a,b,Q,P,ortho)

% Jacobi (for a > -1, b > -1)
% Support: ]-1,1[
% Associated Probability measure: Beta

% Legendre (for a = 0, b = 0)
% Support: ]-1,1[
% Associated Probability measure: Uniform

% Number of experiments
K = size(Q,2); 

% Basis variable preallocation
basis = zeros(P+1,K);

basis(1,:) = 1;                     % P_0(Q)
basis(2,:) = (a-b+(a+b+2)*Q)/2;     % P_1(Q)

for i = 2:P
    A_ell = ((a+b+2*i-1).*((a+b+2*i-2).*(a+b+2*i).*Q + a^2 - b^2)) ./ (2*i.*(a+b+i).*(a+b+2*i-2));
    B_ell = (2*(a+i-1).*(b+i-1).*(a+b+2*i)) ./ (2*i.*(a+b+i).*(a+b+2*i-2));
    basis(i+1,:) = A_ell.*basis(i,:) - B_ell.*basis(i-1,:);   % P_k(Q)
end

if strcmp(ortho,'y')
    if a==0 && b==0     % Legendre
        for i=0:P
            basis(i+1,:) = basis(i+1,:)*sqrt((2*i+1)/2);
        end
    else                % Jacobi
        for i=0:P
            basis(i+1,:) = basis(i+1,:)*...
                sqrt(factorial(i)*(a+b+1+2*i)*gamma(a+b+i+1)/(2^(a+b+1)*gamma(a+i+1)*gamma(b+i+1)));
        end
    end
end