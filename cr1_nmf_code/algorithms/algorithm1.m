% algorithm 1 in the paper
% See Z. Liu and V. Y. F. Tan, "Relative Error Bounds for Nonnegative Matrix
% Factorization under a Geometric Assumption" 
%
% [ Ind ] = algorithm1( V,K )
%
% Input.
%   V              : (F x N) matrix to factorize
%   K              : latent dimensionality
%
% Output.
%   Ind          : final cluster indicator vector of size (N x 1)


function [ Ind ] = algorithm1( V,K )

[F,N] = size(V);
centroids = zeros(F,K);
V_norm = V./repmat(sqrt(sum(V.*V)),F,1);
sel = randi(N,1,1);
centroids(:,1) = V_norm(:,sel);
innerprod_mat = zeros(N,K);
for k = 2:K
    innerprod_mat(:,k-1) = V_norm'*centroids(:,k-1);%slowest part
    maxinnerprod_vec = max(innerprod_mat(:,1:(k-1)),[],2);
    [~,ind_temp] = min(maxinnerprod_vec);
    centroids(:,k) = V_norm(:,ind_temp);
end
innerprod_mat(:,K) = V_norm'*centroids(:,K);
[~,Ind] = max(innerprod_mat,[],2);

