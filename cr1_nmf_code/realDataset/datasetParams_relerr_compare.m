% this function computes the results for Table III: RUNNING TIMES FOR INITIALIZATION
% and Table IV: SHIFT NUMBER FOR INITIALIZATION APPROACHES

% [ iter_mult,iter_hals,shift_mult,shift_hals, running_time ] = datasetParams_time_compare( V,K_est,maxiter,multiTimes )
%
% Input.
%   V              : (F x N) matrix to factorize
%   (W,H)          : initial matrices of dimensions (F x K) and (K x N)
%   K_est          : latent dimensionality
%   maxiter        : maximum number of runs for taking average
%   multiTimes     : the running time of mult and hals is multiTimes times 
%                    the running time of the maximal running time of the 
%                    three initialization methods
%
% Output.
%   iter_mult      : the total number of iterations for mult
%   iter_hals      : the total number of iterations for hals
%   shift_mult     : a vector recording the shifts for the three
%                    initialization methods with respect to mult;
%   shift_hals     : a vector recording the shifts for the three
%                    initialization methods with respect to hals;
%   running_time   : the initialization time for the three initialization 
%                    methods for a dataset

function [ iter_mult,iter_hals,shift_mult,shift_hals, running_time ] = datasetParams_relerr_compare( V,K_est,maxiter,multiTimes )
if nargin < 4
    multiTimes = 30;
    % the running time of mult and hals is multiTimes times the running time of
    % the maximal running time of the three initialization methods
end
if nargin < 3
    maxiter = 20; % to take the average
end

[F,N] = size(V);
running_time = zeros(maxiter,3);
for iter_num = 1:maxiter
    t0 = cputime; [W_cr1nmf,~] = cr1nmf( V,K_est ); 
    running_time(iter_num,1) = cputime-t0;
    
    spkm_maxiter = 10; t0 = cputime; W_spkm = spkm( V,K_est,spkm_maxiter ); 
    running_time(iter_num,2) = cputime-t0;
    
    t0 = cputime; [W_nnd,~] = nndsvd( V,K_est ); 
    running_time(iter_num,3) = cputime-t0;
end
runtime_inits = mean(running_time); % obtain the average running time for the initialization methods

timelimit = floor(max(runtime_inits) * multiTimes) +  1; maxiterAlgs = 1e5; 
iter_mult = 0; iter_hals = 0;
for iter_num = 1:floor(maxiter/4)
    W0 = rand(F,K_est); H0 = rand(K_est,N);
    [~,~,rel_err_mult,~] = mult(V,W0,H0,maxiterAlgs,timelimit); 
    iter_mult = iter_mult+length(rel_err_mult);
    
    [~,~,rel_err_hals,~] = hals(V,W0,H0,maxiterAlgs,timelimit); 
    iter_hals = iter_hals+length(rel_err_hals);
end
iter_mult = floor(iter_mult/floor(maxiter/4));iter_hals=floor(iter_hals/floor(maxiter/4));
runtime_mu = timelimit/iter_mult;runtime_hals=timelimit/iter_hals;

% the running time of initialization methods divided by the running time of
% mult and hals
shift_mult = floor(runtime_inits/runtime_mu); 
shift_hals = floor(runtime_inits/runtime_hals);
