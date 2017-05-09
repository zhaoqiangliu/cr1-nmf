% Demo file for ARD-NMF
% Just run demo.m in Matlab

clear all; clc

%% synthetic experiments to compare cr1-nmf, mult, nnlsb and hals methods
addpath Synthetic/
addpath algorithms/

N = 10000; % number of samples
F = 1600; % dimensionality
K = 40; % number of circular cones
max_alpha = 2*0.2; % 2 times of size angle
pert_angle = 0.01; % perturbation angle
mu_vec = (1:K); % parameters for exponential distributions which generate the squared length of samples
timelimit = inf; % maximal running time
mult_maxiter = 1000; iter_nnlsb = 100; hals_maxiter = 400; % maximal number of iterations for three methods
proj_label = 1; % take projection

[ u_mat,alpha_vec ] = genCones( F,K,max_alpha, pert_angle); % generate the basis vectors and size angles 
while(~isempty(find(u_mat(:)<0, 1)))
    fprintf('There are negative values in basis vectors, run again..\n');
    [ u_mat,alpha_vec ] = genCones( F,K,max_alpha, pert_angle);
end

[ V,~ ] = generate_allpts( u_mat,alpha_vec,mu_vec,N,proj_label ); % generate data matrix
fprintf(['The mean value of length is ',num2str(sum(sum(V.*V))/N),'\n']);

%% run cr1-nmf, mult, nnlsb, hals
normV = norm(V,'fro');

t0 = cputime;
[W_cr1nmf,H_cr1nmf] = cr1nmf( V,K );
rel_err_cr1nmf = norm(W_cr1nmf*H_cr1nmf-V,'fro')/normV; % relative error for cr1-nmf
t_cr1nmf = cputime - t0; % running time for cr1-nmf
fprintf(['The running time of cr1-nmf is ',num2str(t_cr1nmf),'\n']);

W0 = rand(F,K); H0 = rand(K,N);
[~,~,rel_err_mult,t_mult] = mult(V,W0,H0,mult_maxiter,timelimit); % sequence of relative error and running time for mult
ind_mult = find(rel_err_mult <= rel_err_cr1nmf,1);
t_mult = t_mult(ind_mult);
fprintf(['The running time when mult first achieve the same relative error is ',num2str(t_mult),'\n']);

[~,~,rel_err_nnlsb,t_nnlsb] = nnlsb( V,W0,iter_nnlsb,timelimit ); % sequence of relative error and running time for nnlsb
ind_nnlsb = find(rel_err_nnlsb <= rel_err_cr1nmf,1);
t_nnlsb = t_nnlsb(ind_nnlsb);
fprintf(['The running time when nnlsb first achieve the same relative error is ',num2str(t_nnlsb),'\n']);

[~,~,rel_err_hals,t_hals] = hals(V,W0,H0,hals_maxiter,timelimit); % sequence of relative error and running time for hals
ind_hals = find(rel_err_hals <= rel_err_cr1nmf,1);
t_hals = t_hals(ind_hals);
fprintf(['The running time when hals first achieve the same relative error is ',num2str(t_hals),'\n']);


%% real-dataset experiments to see the initialization effect of various methods
% In this toy experiment, for simplicity, we ignore the running time of
% initialization methods when comparing the performance of standard NMF algorithms

%% Load data and set parameters
clear all; clc
X = load('dataset/PaviaU.mat'); V_temp = X.paviaU;
V = reshape(V_temp,size(V_temp,1)*size(V_temp,2),size(V_temp,3)); % data matrix
K_est = 9; % latent dimensionality of NMF
[F,N] = size(V);
normV = norm(V,'fro');

% maximal number of iterations of standard NMF algorithms
mult_maxiter = 800; 
nnlsb_maxiter = 50;
hals_maxiter = 200;

spkm_maxiter = 10; % number of iterations for spkm initialization
timelimit = 1e6;

%% initializations
t0 = cputime; 
[W_cr1nmf,H_cr1nmf] = cr1nmf( V,K_est ); 
epsilon_cr1nmf = mean(H_cr1nmf(H_cr1nmf>0));
t_cr1nmf = cputime-t0; % cr1-nmf initialization
rel_err_cr1nmf = norm(W_cr1nmf*H_cr1nmf-V,'fro')/normV;

t0 = cputime; 
W_spkm = spkm( V,K_est,spkm_maxiter ); 
H_spkm = rand(K_est,N);
t_spkm = cputime-t0; % spkm initialization

t0 = cputime; 
[W_nndsvd,H_nndsvd] = nndsvd( V,K_est );
t_nndsvd = cputime-t0; % nndsvd initialization

fprintf(['The running time of cr1-nmf is ',num2str(t_cr1nmf),'\n']);
fprintf(['The running time of spkm is ',num2str(t_spkm),'\n']);
fprintf(['The running time of nndsvd is ',num2str(t_nndsvd),'\n']);

fprintf(['The relative error of cr1-nmf is ',num2str(rel_err_cr1nmf),'\n']);

%% combining with standard NMF algorithms
W0 = rand(F,K_est); H0 = rand(K_est,N);
   
% combining with mult
[~,~,rel_err_mult,t_mult] = mult(V,W0,H0,mult_maxiter,timelimit); % random initialization
[~,~,rel_err_cr1nmf_mult,t_cr1nmf_mult] = mult(V,W_cr1nmf,H_cr1nmf+epsilon_cr1nmf*ones(size(H_cr1nmf)),...
    mult_maxiter,timelimit); % cr1-nmf initialization
[~,~,rel_err_spkm_mult,t_spkm_mult] = mult(V,W_spkm,H_spkm,mult_maxiter,timelimit); % spkm initialization
[~,~,rel_err_nndsvd_mult,t_nndsvd_mult] = mult(V,W_nndsvd,H_nndsvd,mult_maxiter,timelimit); % nndsvd initialization

% combining with nnlsb
[~,~,rel_err_nnlsb,t_nnlsb] = nnlsb( V,W0,nnlsb_maxiter,timelimit ); % random initialization
[~,~,rel_err_cr1nmf_nnlsb,t_cr1nmf_nnlsb] = nnlsb( V,W_cr1nmf,nnlsb_maxiter,timelimit ); % cr1-nmf initialization
[~,~,rel_err_spkm_nnlsb,t_spkm_nnlsb] = nnlsb( V,W_spkm,nnlsb_maxiter,timelimit ); % spkm initialization
[~,~,rel_err_nndsvd_nnlsb,t_nndsvd_nnlsb] = nnlsb( V,W_nndsvd,nnlsb_maxiter,timelimit ); % nndsvd initialization

% combining with hals
[~,~,rel_err_hals,t_hals] = hals(V,W0,H0,hals_maxiter,timelimit); % random initialization
[~,~,rel_err_cr1nmf_hals,t_cr1nmf_hals] = hals(V,W_cr1nmf,H_cr1nmf+epsilon_cr1nmf*ones(size(H_cr1nmf)),...
    hals_maxiter,timelimit); % cr1-nmf initialization
[~,~,rel_err_spkm_hals,t_spkm_hals] = hals(V,W_spkm,H_spkm,hals_maxiter,timelimit); % spkm initialization
[~,~,rel_err_nndsvd_hals,t_nndsvd_hals] = hals(V,W_nndsvd,H_nndsvd,hals_maxiter,timelimit); % nndsvd initialization

%% Display the results
addpath plotFuncs/
str_now=datestr(now,30);
dataset_name = 'PaviaU';

% mult
h1 = figure('position',[0 0 300 360]);
plot(rel_err_mult,'r--', 'Linewidth', 1.5);hold on;
plot(rel_err_cr1nmf_mult,'b', 'Linewidth', 1.5);hold on;
plot(rel_err_spkm_mult,'g', 'Linewidth', 1.5);hold on;
plot(rel_err_nndsvd_mult,'k-.', 'Linewidth', 1.5);hold on;
xlabel('Iteration Number','FontSize',18);
ylabel('Relative Error','FontSize',18);
title('MULT', 'Fontsize',18);
legend('rand','cr1-nmf','spkm','nndsvd','Location', 'Northeast');
set(gcf,'color','w');set(gca,'FontSize',18);
filename=['Figures/realDataset/MULT_dataset_',dataset_name,'_K_est',num2str(K_est),'_',str_now(1:8),'.pdf'];
export_fig(gcf,'Color','Transparent',filename);

% nnlsb
h2 = figure('position',[0 0 300 360]);
plot(rel_err_nnlsb,'r--', 'Linewidth', 1.5);hold on;
plot(rel_err_cr1nmf_nnlsb,'b', 'Linewidth', 1.5);hold on;
plot(rel_err_spkm_nnlsb,'g', 'Linewidth', 1.5);hold on;
plot(rel_err_nndsvd_nnlsb,'k-.', 'Linewidth', 1.5);hold on;
xlabel('Iteration Number','FontSize',18);
ylabel('Relative Error','FontSize',18);
title('NNLSB', 'Fontsize',18);
legend('rand','cr1-nmf','spkm','nndsvd','Location', 'Northeast');
set(gcf,'color','w');set(gca,'FontSize',18);
filename=['Figures/realDataset/NNLSB_dataset_',dataset_name,'_K_est',num2str(K_est),'_',str_now(1:8),'.pdf'];
export_fig(gcf,'Color','Transparent',filename);

% hals
h3 = figure('position',[0 0 300 360]);
plot(rel_err_hals,'r--', 'Linewidth', 1.5);hold on;
plot(rel_err_cr1nmf_hals,'b', 'Linewidth', 1.5);hold on;
plot(rel_err_spkm_hals,'g', 'Linewidth', 1.5);hold on;
plot(rel_err_nndsvd_hals,'k-.', 'Linewidth', 1.5);hold on;
xlabel('Iteration Number','FontSize',18);
ylabel('Relative Error','FontSize',18);
title('HALS', 'Fontsize',18);
legend('rand','cr1-nmf','spkm','nndsvd','Location', 'Northeast');
set(gcf,'color','w');set(gca,'FontSize',18);
filename=['Figures/realDataset/HALS_dataset_',dataset_name,'_K_est',num2str(K_est),'_',str_now(1:8),'.pdf'];
export_fig(gcf,'Color','Transparent',filename);
