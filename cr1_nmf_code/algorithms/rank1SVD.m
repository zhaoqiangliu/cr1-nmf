% This function is used in cr1nmf.m
function [ u,s,v ] = rank1SVD( A,maxit,tol )
% calculate the rank one SVD of matrix A
if nargin<3
    tol=1e-4;
end
if nargin<2
    maxit=100;
end

[F,N]=size(A);
if F/N>1.0713
    ddata = A'*A;
    ddata = max(ddata,ddata');
    [ v, s_sq, ~ ] = power_method (  ddata, rand(N,1), maxit, tol );
    clear ddata;
    s=sqrt(s_sq);
    u=A*v/s;
else
    ddata = A*A';
    ddata = max(ddata,ddata');
    [ u, s_sq, ~ ] = power_method (  ddata, rand(F,1), maxit, tol );
    s=sqrt(s_sq);
    clear ddata;
    if(nargout>2)
        v=A'*u/s;
    end
end





