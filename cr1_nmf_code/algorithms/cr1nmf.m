% cr1-nmf algorithm
% See Z. Liu and V. Y. F. Tan, "Relative Error Bounds for Nonnegative Matrix
% Factorization under a Geometric Assumption" 
%
% [W,H]  = cr1nmf( V,K )
%
% Input.
%   V              : (F x N) matrix to factorize
%   K              : latent dimensionality
%
% Output.
%   (W,H)          : final nonnegative matrices of dimensions (F x K) and (K x N)


function [W,H]  = cr1nmf( V,K )
[F,N] = size(V);
[ Ind ] = algorithm1( V,K );

W = zeros(F,K);
H = zeros(K,N);
for k=1:K
    Vk=V(:,Ind == k);
    [u,s,v] = rank1SVD(Vk);
    W(:,k) = abs(u);hk = s*abs(v);
    H(k,Ind == k) = hk;
end







