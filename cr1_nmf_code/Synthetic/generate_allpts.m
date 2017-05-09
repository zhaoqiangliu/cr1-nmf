% this function is used to generate N_all points according to parameters
%
% [ V_all,cone_inds,P_mat ] = generate_allpts( u_mat,alpha_vec,mu_vec,N_all,proj_label )
%
% Input.
%   u_mat          : (F x K) matrix; contains the basis vectors;
%   alpha_vec      : (1 x K) vector; contains the size angles;
%   mu_vec         : (1 x K) vector; mu_vec = 1./lambda_vec, where lambda_vec 
%                    is the parameters for exponential distributions (cf. Step 2 in the generative process)
%   N_all          : total number of data points to generate
%   proj_label     : proj_label = 1 means taking projection to nonnegative orthant (Step 4 in the generative process)
%
% Output.
%   V_all          : (F x N_all) matrix; the data matrix generated
%   cone_inds      : (1 x N_all) vector; containing values in {1,2,...,K}, recording which
%                    cone the data points were generated from
%   P_mat          : (KF x F) matrix; recording the K Householder transform matrices
%                    corresponding to the K basis vectors

function [ V_all,cone_inds,P_mat ] = generate_allpts( u_mat,alpha_vec,mu_vec,N_all,proj_label )

[F,K]=size(u_mat);
cone_inds=randi(K,1,N_all);
num_vec=zeros(1,K);
for k=1:K
    num_vec(k)=sum(cone_inds==k);
end

len_mat=zeros(K,max(num_vec));
for k=1:K
    len_mat(k,1:num_vec(k))=exprnd(mu_vec(k),1,num_vec(k));
end
len_mat=sqrt(len_mat);

angle_mat=zeros(K,max(num_vec));
for k=1:K
    angle_mat(k,1:num_vec(k))=unifrnd(0,alpha_vec(k),1,num_vec(k));
end

val_mat=normrnd(0,1,F-1,N_all);
val_col_sqrtSum=sqrt(sum(val_mat.*val_mat));
val_col_sqrtSumMat=repmat(val_col_sqrtSum,F-1,1);
null_mat=val_mat./val_col_sqrtSumMat;

Kcounts=zeros(1,K);V_all=zeros(F,N_all);
e1=zeros(F,1);e1(1)=1;
u_mat_origin=repmat(e1,1,K);
w_mat=u_mat_origin-u_mat;
w_mat=w_mat./repmat(sqrt(sum(w_mat.*w_mat)),F,1);

for n=1:N_all
    nk=cone_inds(n);Kcounts(nk)=Kcounts(nk)+1;
    V_all(:,n)=len_mat(nk,Kcounts(nk))*(u_mat_origin(:,nk)*cos(angle_mat(nk,Kcounts(nk)))+sin(angle_mat(nk,Kcounts(nk)))*[0;null_mat(1:end,n)]);
end

for k=1:K
    P_mat(((k-1)*F+1):(k*F),:)=eye(F)-2*w_mat(:,k)*w_mat(:,k)';
    V_all(:,cone_inds==k)=P_mat(((k-1)*F+1):(k*F),:)*V_all(:,cone_inds==k);
end

% projection to the nonnegative orthant
if(proj_label==1)
    V_all_len=sqrt(sum(V_all.*V_all));
    for n=1:N_all
        if(sum(V_all(:,n)<0)>0)
            V_all(V_all(:,n)<0,n)=0;
            V_all(:,n)=V_all_len(n)*V_all(:,n)/norm(V_all(:,n));
        end
    end
end