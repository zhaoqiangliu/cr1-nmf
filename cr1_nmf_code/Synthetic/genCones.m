function [ u_mat,alpha_vec ] = genCones( F,K,max_alpha, pert_angle)
% this function is used to generate the matrix of basis vectors, u_mat;
% and the vector of size angles, alpha_vec

alpha_vec=0.5*max_alpha*ones(1,K);
alpha_mat=zeros(K);
for i=1:K
    for j=(i+1):K
        alpha_mat(i,j)=max(3*alpha_vec(i)+alpha_vec(j),3*alpha_vec(j)+alpha_vec(i))+pert_angle;
        alpha_mat(j,i)=alpha_mat(i,j);
    end
end
u1=rand(F,1);u1=u1+0.5;u1=u1/norm(u1);
u_mat=zeros(F,K);u_mat(:,1)=u1;

for k=2:K
    A=u_mat(:,1:(k-1))';
    y=cos(alpha_mat(1:(k-1),k));
    u_ln=A'*((A*A')\y);len=sqrt(1-norm(u_ln)^2);
    pos=find(u_ln==min(u_ln));
    AA=A(:,setdiff(1:F,pos));
    yy=-A(:,pos);u_lnln=AA'*((AA*AA')\yy);
    u_out=[u_lnln(1:(pos-1));1;u_lnln(pos:end)];
    u_out=len*u_out/norm(u_out);
    u=u_ln+u_out;
    u_mat(:,k)=u;
end

