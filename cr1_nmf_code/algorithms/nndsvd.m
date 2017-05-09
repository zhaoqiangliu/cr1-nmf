function [ W,H ] = nndsvd( V,K )
[F,N]=size(V);

% we use the accelerated singular value decomposition
[U, S, VV] = mySVD(V,K); 

W=zeros(F,K);H=zeros(K,N);
W(:,1)=sqrt(S(1,1))*abs(U(:,1));
H(1,:)=sqrt(S(1,1))*abs(VV(:,1)');

for k=2:K
    x=U(:,k);y=VV(:,k);
    xp=x.*(x>=0);xn=(-x).*(x<0);
    yp=y.*(y>=0);yn=(-y).*(y<0);
    
    xpnrm=norm(xp);ypnrm=norm(yp);mp=xpnrm*ypnrm;
    xnnrm=norm(xn);ynnrm=norm(yn);mn=xnnrm*ynnrm;
    
    if(mp>mn)
        u=xp/xpnrm;v=yp/ypnrm;sigma=mp;
    else
        u=xn/xnnrm;v=yn/ynnrm;sigma=mn;
    end
    
    W(:,k)=sqrt(S(k, k)*sigma)*u;
    H(k,:)=sqrt(S(k, k)*sigma)*v';
end

