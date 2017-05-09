% Initialization Performance in Terms of Clustering

function initialization_clustering_compare(dataset_name)

addpath ../algorithms/;
addpath ../plotFuncs/;

str_now = datestr(now,30);
maxiter = 10; % total number of runs
default_maxiter = 1000; % maximal iteration number of k-means and spkm

% the NMF algorithms are terminated if the variation of the product of factor matrices is small over 10 iterations
early_stopping = 1; 
fprintf(['Now running for dataset ',dataset_name,'...\n']);
[ V,K_est,cluster_gd ] = readDataset_clustering_compare( dataset_name );
fprintf('Dataset obtained.\n');
[F,N] = size(V);
normV = norm(V,'fro');

mult_maxiter = 2 * default_maxiter; % because of mult's slow convergence rate
hals_maxiter = default_maxiter; % the maximal number of iteration of hals
nnlsb_maxiter = 100; % because of nnlsb's fast convergence rate
rel_err_mult_mat = zeros(4,mult_maxiter);
rel_err_hals_mat = zeros(4,hals_maxiter);
rel_err_nnlsb_mat = zeros(4,nnlsb_maxiter);

spkm_maxiter_init=10; % the maximal number of iteration of spkm when it is used as an initialization method
spkm_maxiter = default_maxiter; % the maximal number of iteration of spkm
kmeans_maxiter = default_maxiter; % the maximal number of iteration of k-means

timelimit = 1e6;
% save clustering performance results for the 14 methods
nmi_vec = zeros(maxiter,14); dice_vec = zeros(maxiter,14); purity_vec = zeros(maxiter,14); 
rel_err_cr1nmf = zeros(1,maxiter); 

for iter_num=1:maxiter
    % initializations
    fprintf(['The ',num2str(iter_num),'th iteration..\n']);
    [W_cr1nmf,H_cr1nmf] = cr1nmf( V,K_est ); epsilon_cr1nmf = mean(H_cr1nmf(H_cr1nmf>0)); % cr1-nmf initialization
    rel_err_cr1nmf(iter_num)=norm(W_cr1nmf*H_cr1nmf-V,'fro')/normV;
    W_spkm = spkm( V,K_est,spkm_maxiter_init ); H_spkm = rand(K_est,N); % spkm initialization
    [W_nndsvd,H_nndsvd] = nndsvd( V,K_est ); % nndsvd initialization (no need to average over different runs);
    
    %% relative error results
    W0 = rand(F,K_est); H0 = rand(K_est,N);
    % combining with mult algorithm
    [W_mult,H_mult,rel_err_mult] = mult(V,W0,H0,mult_maxiter,timelimit,early_stopping); 
    rel_err_mult_mat(1,1:length(rel_err_mult)) = rel_err_mult_mat(1,1:length(rel_err_mult))+rel_err_mult; 
    rel_err_mult_mat(1,1:length(rel_err_mult)) = 0.5 * rel_err_mult_mat(1,1:length(rel_err_mult)); % random initialization
    
    [W_cr1nmf_mult,H_cr1nmf_mult,rel_err_cr1nmf_mult] = mult(V,W_cr1nmf,H_cr1nmf+epsilon_cr1nmf*ones(size(H_cr1nmf)),mult_maxiter,timelimit,early_stopping);
    rel_err_mult_mat(2,1:length(rel_err_cr1nmf_mult)) = rel_err_mult_mat(2,1:length(rel_err_cr1nmf_mult))+rel_err_cr1nmf_mult; 
    rel_err_mult_mat(2,1:length(rel_err_cr1nmf_mult)) = 0.5*rel_err_mult_mat(2,1:length(rel_err_cr1nmf_mult)); % cr1-nmf initialization
    
    [W_spkm_mult,H_spkm_mult,rel_err_spkm_mult] = mult(V,W_spkm,H_spkm,mult_maxiter,timelimit,early_stopping);
    rel_err_mult_mat(3,1:length(rel_err_spkm_mult)) = rel_err_mult_mat(3,1:length(rel_err_spkm_mult))+rel_err_spkm_mult; 
    rel_err_mult_mat(3,1:length(rel_err_spkm_mult)) = 0.5*rel_err_mult_mat(3,1:length(rel_err_spkm_mult)); % spkm initialization
    
    [W_nndsvd_mult,H_nndsvd_mult,rel_err_nndsvd_mult] = mult(V,W_nndsvd,H_nndsvd,mult_maxiter,timelimit,early_stopping);
    rel_err_mult_mat(4,1:length(rel_err_nndsvd_mult)) = rel_err_mult_mat(4,1:length(rel_err_nndsvd_mult))+rel_err_nndsvd_mult; 
    rel_err_mult_mat(4,1:length(rel_err_nndsvd_mult)) = 0.5*rel_err_mult_mat(4,1:length(rel_err_nndsvd_mult)); % nndsvd initialization
    
    % combining with nnlsb algorithm
    [W_nnlsb,H_nnlsb,rel_err_nnlsb] = nnlsb( V,W0,nnlsb_maxiter,timelimit );
    rel_err_nnlsb_mat(1,1:length(rel_err_nnlsb)) = rel_err_nnlsb_mat(1,1:length(rel_err_nnlsb))+rel_err_nnlsb; 
    rel_err_nnlsb_mat(1,1:length(rel_err_nnlsb)) = 0.5*rel_err_nnlsb_mat(1,1:length(rel_err_nnlsb)); % random initialization
    
    [W_cr1nmf_nnlsb,H_cr1nmf_nnlsb,rel_err_cr1nmf_nnlsb] = nnlsb( V,W_cr1nmf,nnlsb_maxiter,timelimit,early_stopping );
    rel_err_nnlsb_mat(2,1:length(rel_err_cr1nmf_nnlsb)) = rel_err_nnlsb_mat(2,1:length(rel_err_cr1nmf_nnlsb))+rel_err_cr1nmf_nnlsb;
    rel_err_nnlsb_mat(2,1:length(rel_err_cr1nmf_nnlsb)) = 0.5*rel_err_nnlsb_mat(2,1:length(rel_err_cr1nmf_nnlsb)); % cr1-nmf initialization
    
    [W_spkm_nnlsb,H_spkm_nnlsb,rel_err_spkm_nnlsb] = nnlsb( V,W_spkm,nnlsb_maxiter,timelimit,early_stopping );
    rel_err_nnlsb_mat(3,1:length(rel_err_spkm_nnlsb)) = rel_err_nnlsb_mat(3,1:length(rel_err_spkm_nnlsb))+rel_err_spkm_nnlsb;
    rel_err_nnlsb_mat(3,1:length(rel_err_spkm_nnlsb)) = 0.5*rel_err_nnlsb_mat(3,1:length(rel_err_spkm_nnlsb)); % spkm initialization
    
    [W_nndsvd_nnlsb,H_nndsvd_nnlsb,rel_err_nndsvd_nnlsb] = nnlsb( V,W_nndsvd,nnlsb_maxiter,timelimit,early_stopping );
    rel_err_nnlsb_mat(4,1:length(rel_err_nndsvd_nnlsb)) = rel_err_nnlsb_mat(4,1:length(rel_err_nndsvd_nnlsb))+rel_err_nndsvd_nnlsb;
    rel_err_nnlsb_mat(4,1:length(rel_err_nndsvd_nnlsb)) = 0.5*rel_err_nnlsb_mat(4,1:length(rel_err_nndsvd_nnlsb)); % nndsvd initialization
    
    % combining with hals algorithm
    [W_hals,H_hals,rel_err_hals] = hals(V,W0,H0,hals_maxiter,timelimit,early_stopping);
    rel_err_hals_mat(1,1:length(rel_err_hals)) = rel_err_hals_mat(1,1:length(rel_err_hals))+rel_err_hals;
    rel_err_hals_mat(1,1:length(rel_err_hals)) = 0.5*rel_err_hals_mat(1,1:length(rel_err_hals)); % random initialization
    
    [W_cr1nmf_hals,H_cr1nmf_hals,rel_err_cr1nmf_hals] = hals(V,W_cr1nmf,H_cr1nmf+epsilon_cr1nmf*ones(size(H_cr1nmf)),hals_maxiter,timelimit,early_stopping);
    rel_err_hals_mat(2,1:length(rel_err_cr1nmf_hals)) = rel_err_hals_mat(2,1:length(rel_err_cr1nmf_hals))+rel_err_cr1nmf_hals;
    rel_err_hals_mat(2,1:length(rel_err_cr1nmf_hals)) = 0.5*rel_err_hals_mat(2,1:length(rel_err_cr1nmf_hals)); % cr1-nmf initializations
    
    [W_spkm_hals,H_spkm_hals,rel_err_spkm_hals] = hals(V,W_spkm,H_spkm,hals_maxiter,timelimit,early_stopping);
    rel_err_hals_mat(3,1:length(rel_err_spkm_hals)) = rel_err_hals_mat(3,1:length(rel_err_spkm_hals))+rel_err_spkm_hals;
    rel_err_hals_mat(3,1:length(rel_err_spkm_hals)) = 0.5*rel_err_hals_mat(3,1:length(rel_err_spkm_hals)); % spkm initialization
    
    [W_nndsvd_hals,H_nndsvd_hals,rel_err_nndsvd_hals] = hals(V,W_nndsvd,H_nndsvd,hals_maxiter,timelimit,early_stopping);
    rel_err_hals_mat(4,1:length(rel_err_nndsvd_hals)) = rel_err_hals_mat(4,1:length(rel_err_nndsvd_hals))+rel_err_nndsvd_hals;
    rel_err_hals_mat(4,1:length(rel_err_nndsvd_hals)) = 0.5*rel_err_hals_mat(4,1:length(rel_err_nndsvd_hals)); % nndsvd initialization
    
    
    % clustering results
    [ind_kmeans,~] = litekmeans(V',K_est,'MaxIter',kmeans_maxiter);
    nmi_vec(iter_num,1) = nmi(ind_kmeans',cluster_gd);
    dice_vec(iter_num,1) = Dice(ind_kmeans',cluster_gd); 
    purity_vec(iter_num,1) = purity(ind_kmeans',cluster_gd); % k-means, we use a faster version of k-means
    
    [~,ind_spkm] = spkm( V,K_est,spkm_maxiter );
    nmi_vec(iter_num,2) = nmi(ind_spkm,cluster_gd);
    dice_vec(iter_num,2) = Dice(ind_spkm,cluster_gd); 
    purity_vec(iter_num,2) = purity(ind_spkm,cluster_gd); % spkm
    
    % clustering performance for mult algorithms
    norm_vec = sqrt(sum(W_mult.*W_mult)); H_mult = diag(norm_vec)*H_mult;
    [~,ind_mult]=max(H_mult);
    nmi_vec(iter_num,3) = nmi(ind_mult,cluster_gd);
    dice_vec(iter_num,3) = Dice(ind_mult,cluster_gd); 
    purity_vec(iter_num,3) = purity(ind_mult,cluster_gd); % mult
    
    norm_vec = sqrt(sum(W_cr1nmf_mult.*W_cr1nmf_mult)); H_cr1nmf_mult = diag(norm_vec)*H_cr1nmf_mult;
    [~,ind_cr1nmf_mult] = max(H_cr1nmf_mult);
    nmi_vec(iter_num,4) = nmi(ind_cr1nmf_mult,cluster_gd);
    dice_vec(iter_num,4) = Dice(ind_cr1nmf_mult,cluster_gd); 
    purity_vec(iter_num,4) = purity(ind_cr1nmf_mult,cluster_gd); % cr1-nmf + mult
    
    norm_vec = sqrt(sum(W_spkm_mult.*W_spkm_mult)); H_spkm_mult = diag(norm_vec)*H_spkm_mult;
    [~,ind_spkm_mult] = max(H_spkm_mult);
    nmi_vec(iter_num,5) = nmi(ind_spkm_mult,cluster_gd);
    dice_vec(iter_num,5) = Dice(ind_spkm_mult,cluster_gd);
    purity_vec(iter_num,5) = purity(ind_spkm_mult,cluster_gd); % spkm + mult
    
    norm_vec = sqrt(sum(W_nndsvd_mult.*W_nndsvd_mult)); H_nndsvd_mult = diag(norm_vec)*H_nndsvd_mult;
    [~,ind_nndsvd_mult] = max(H_nndsvd_mult);
    nmi_vec(iter_num,6) = nmi(ind_nndsvd_mult,cluster_gd);
    dice_vec(iter_num,6) = Dice(ind_nndsvd_mult,cluster_gd);
    purity_vec(iter_num,6)=purity(ind_nndsvd_mult,cluster_gd); % nndsvd + mult
    
    % clustering performance for nnlsb algorithms
    norm_vec = sqrt(sum(W_nnlsb.*W_nnlsb)); H_nnlsb = diag(norm_vec)*H_nnlsb;
    [~,ind_nnlsb] = max(H_nnlsb);
    nmi_vec(iter_num,7) = nmi(ind_nnlsb,cluster_gd);
    dice_vec(iter_num,7) = Dice(ind_nnlsb,cluster_gd);
    purity_vec(iter_num,7) = purity(ind_nnlsb,cluster_gd); % nnlsb
    
    norm_vec = sqrt(sum(W_cr1nmf_nnlsb.*W_cr1nmf_nnlsb)); H_cr1nmf_nnlsb = diag(norm_vec)*H_cr1nmf_nnlsb;
    [~,ind_cr1nmf_nnlsb] = max(H_cr1nmf_nnlsb);
    nmi_vec(iter_num,8) = nmi(ind_cr1nmf_nnlsb,cluster_gd);
    dice_vec(iter_num,8) = Dice(ind_cr1nmf_nnlsb,cluster_gd);
    purity_vec(iter_num,8) = purity(ind_cr1nmf_nnlsb,cluster_gd); % cr1-nmf + nnlsb
    
    norm_vec = sqrt(sum(W_spkm_nnlsb.*W_spkm_nnlsb)); H_spkm_nnlsb = diag(norm_vec)*H_spkm_nnlsb;
    [~,ind_spkm_nnlsb] = max(H_spkm_nnlsb);
    nmi_vec(iter_num,9) = nmi(ind_spkm_nnlsb,cluster_gd);
    dice_vec(iter_num,9) = Dice(ind_spkm_nnlsb,cluster_gd);
    purity_vec(iter_num,9) = purity(ind_spkm_nnlsb,cluster_gd); % spkm + nnlsb
    
    norm_vec = sqrt(sum(W_nndsvd_nnlsb.*W_nndsvd_nnlsb)); H_nndsvd_nnlsb = diag(norm_vec)*H_nndsvd_nnlsb;
    [~,ind_nndsvd_nnlsb] = max(H_nndsvd_nnlsb);
    nmi_vec(iter_num,10) = nmi(ind_nndsvd_nnlsb,cluster_gd);
    dice_vec(iter_num,10) = Dice(ind_nndsvd_nnlsb,cluster_gd);
    purity_vec(iter_num,10) = purity(ind_nndsvd_nnlsb,cluster_gd); % nndsvd + nnlsb
    
    % clustering performance for hals algorithms
    norm_vec = sqrt(sum(W_hals.*W_hals)); H_hals = diag(norm_vec)*H_hals;
    [~,ind_hals] = max(H_hals);
    nmi_vec(iter_num,11) = nmi(ind_hals,cluster_gd);
    dice_vec(iter_num,11) = Dice(ind_hals,cluster_gd);
    purity_vec(iter_num,11) = purity(ind_hals,cluster_gd); % hals
    
    norm_vec = sqrt(sum(W_cr1nmf_hals.*W_cr1nmf_hals)); H_cr1nmf_hals = diag(norm_vec)*H_cr1nmf_hals;
    [~,ind_cr1nmf_hals] = max(H_cr1nmf_hals);
    nmi_vec(iter_num,12) = nmi(ind_cr1nmf_hals,cluster_gd);
    dice_vec(iter_num,12) = Dice(ind_cr1nmf_hals,cluster_gd);
    purity_vec(iter_num,12) = purity(ind_cr1nmf_hals,cluster_gd); % cr1-nmf + hals
    
    norm_vec = sqrt(sum(W_spkm_hals.*W_spkm_hals)); H_spkm_hals = diag(norm_vec)*H_spkm_hals;
    [~,ind_spkm_hals] = max(H_spkm_hals);
    nmi_vec(iter_num,13) = nmi(ind_spkm_hals,cluster_gd);
    dice_vec(iter_num,13) = Dice(ind_spkm_hals,cluster_gd);
    purity_vec(iter_num,13) = purity(ind_spkm_hals,cluster_gd); % spkm + hals
    
    norm_vec = sqrt(sum(W_nndsvd_hals.*W_nndsvd_hals)); H_nndsvd_hals = diag(norm_vec)*H_nndsvd_hals;
    [~,ind_nndsvd_hals] = max(H_nndsvd_hals);
    nmi_vec(iter_num,14) = nmi(ind_nndsvd_hals,cluster_gd);
    dice_vec(iter_num,14) = Dice(ind_nndsvd_hals,cluster_gd);
    purity_vec(iter_num,14) = purity(ind_nndsvd_hals,cluster_gd); % nndsvd + hals
end

filename = ['../output/realDataset/rel_err_mult_mat_dataset_',dataset_name,'_K_est',num2str(K_est),'_',str_now(1:8),'.mat'];
save(filename,'rel_err_mult_mat');
filename = ['../output/realDataset/rel_err_nnlsb_mat_dataset_',dataset_name,'_K_est',num2str(K_est),'_',str_now(1:8),'.mat'];
save(filename,'rel_err_nnlsb_mat');
filename = ['../output/realDataset/rel_err_hals_mat_dataset_',dataset_name,'_K_est',num2str(K_est),'_',str_now(1:8),'.mat'];
save(filename,'rel_err_hals_mat');
nmidicepurity_mat = [nmi_vec;dice_vec;purity_vec];
filename = ['../output/realDataset/nmidicepurity_mat_dataset_',dataset_name,'_K_est',num2str(K_est),'_',str_now(1:8),'.mat'];
save(filename,'nmidicepurity_mat');
filename = ['../output/realDataset/rel_err_cr1nmf_dataset_',dataset_name,'_K_est',num2str(K_est),'_',str_now(1:8),'.mat'];
save(filename,'rel_err_cr1nmf');