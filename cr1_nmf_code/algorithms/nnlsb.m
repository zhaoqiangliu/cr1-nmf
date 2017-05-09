% the alternating nonnegative least-squares algorithm with block-pivoting
% written according to the paper: J. Kim and H. Park. Toward faster 
% nonnegative matrix factorization: A new algorithm and comparisons.
% 
% [ W,H,err,t,iter ] = nnlsb( V,W,maxiter,timelimit,early_stopping,tol )
%
% Input.
%   V              : (F x N) matrix to factorize
%   W              : initial matrix of dimension (F x K); note that for
%                    this algorithm, we only need to initialize W (or H)
%   maxiter        : maximum number of iterations
%   timelimit      : maximum time alloted to the algorithm
%   early_stopping : if early_stopping == 1, terminate the algorithm when the
%                    variation of product of factor matrices (<tol) is small 
%                    over 10 times
%   tol            : used for early_stopping
%
% Output.
%   (W,H)          : final nonnegative matrices, WH approximate V
%   (rel_err,t)    : relative error and time after each iteration
%   iter           : number of iterations when terminated

function [ W,H,rel_err,t,iter ] = nnlsb( V,W,maxiter,timelimit,early_stopping,tol )

if nargin <= 2, maxiter = 20; end
if nargin <= 3, timelimit = inf; end
if nargin <= 4, early_stopping = 0; end 
if nargin <= 5, tol = 1e-4; end 

etime = cputime;
normV = norm(V,'fro'); nV = normV*normV;
err = []; t = []; iter = 0;
while iter < maxiter && cputime - etime <= timelimit
    H = nnlsm_blockpivot( W,V,0);
    W_trans = nnlsm_blockpivot( H',V',0);
    W = W_trans';
    
    A = (W'*V); B = (W'*W); 
    cnT = cputime;
    err = [err sqrt( (nV-2*sum(sum(H.*A)) + sum(sum(B.*(H*H')))) )];
    etime = etime+(cputime-cnT);
    t = [t cputime-etime];
    if(length(err) > 9)
        temp_e = err((end-9):end); temp_err = abs(temp_e(1)-temp_e(end));
        if(early_stopping == 1 && temp_err < tol)
            fprintf(['nnlsb algorithm converged, total iteration number is ', num2str(iter+1),'\n']);
            break;
        end
    end
    iter = iter+1;
end
rel_err = err/normV;


