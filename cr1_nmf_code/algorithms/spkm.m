% spherical k-means algorithm
%
% [ u_init, Y_old,iter ] = spkm( V,K,maxiter,norm_label )
%
% Input.
%   V              : (F x N) matrix to factorize
%   K              : number of clusters
%   maxiter        : maximum number of iterations
%   norm_label     : if norm_label == 0, we need to normalize V first
%
% Output.
%   u_init         : (F x K) matrix; containing the centroids
%   Y_old          : (1 x N) vector; containing the cluster indicator
%                    values for each data point
%   iter           : number of iterations when terminated

function [ u_init, Y_old,iter ] = spkm( V,K,maxiter,norm_label )

[F,N]=size(V);
if nargin<= 3, norm_label=0; end
if(norm_label==0), V=V./repmat(sqrt(sum(V.*V)),F,1); end

% initialization of centroids
Krand=randperm(N);Krand=Krand(1:K);
u_init=V(:,Krand);

Y_old=zeros(1,N);
for iter=1:maxiter
    temp_mat=u_init'*V;
    [~,Y_new]=max(temp_mat);
    count=K-length(unique(Y_new));null_row=setdiff(1:K,unique(Y_new));
    
    % To deal with empty clusters. We use the strategy given in S. Zhong, 
    % Efficient Online Spherical K-means Clustering
    % That is, At the end of each batch iteration (i.e., each time after
    % going through all data points), we keep a list of up to K data
    % points that are most distant from their cluster centroids, and
    % a list of empty clusters. The centroid of each empty cluster
    % will then be replaced by a data point (from the list constructed
    % above) that is most remote to its center
    if(count>0)
        for k=1:K
            temp_mat(k,Y_new~=k)=0;
        end
        ind_pos=find(temp_mat(:)>0);
        vec_pos=temp_mat(temp_mat>0);
        index_count_pos=smallestK( vec_pos,count );
        index_count_pos=ind_pos(index_count_pos);
        for count_num=1:count
            [~,temp_N]=ind2sub(size(temp_mat),index_count_pos(count_num));
            Y_new(temp_N)=null_row(count_num);
        end
    end
    
    E = sparse(1:N,Y_new,1,N,K,N);
    u_init = V*E;
    u_init = u_init.*repmat((sum(u_init.^2)+1e-16).^(-0.5),F,1);
    
    if(sum(Y_old==Y_new)==N)
        % if the clustering membership does not change, terminate the program
        % fprintf(['spkm algorithm converged, total iteration number is ', num2str(iter),'\n'])
        break;
    else
        Y_old=Y_new;
    end
end

% find the indices for the K smallest entries
function [ index_vec ] = smallestK( vec,K )
index_vec=zeros(1,K);
for k=1:K
    [~,index_vec(k)]=min(vec);
    vec(index_vec(k))=inf;
end
