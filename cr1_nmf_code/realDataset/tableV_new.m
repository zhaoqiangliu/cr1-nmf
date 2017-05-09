% This function produce the figures to illustrate Table V into graphs

function tableV_new(dataset_name)

addpath ../algorithms/;
addpath ../plotFuncs/;

spkm_maxiter = 10;
nnlsb_maxiter = inf;
str_now=datestr(now,30);


fprintf(['Now running for dataset ',dataset_name,'...\n']);
% obtaining parameters for the dataset
[ V,K_est,timelimit ]=readDataset_nnlsb(dataset_name);
fprintf('The dataset is obtained.\n');
[F,~]=size(V);

%initializations
t0 = cputime; [W_cr1nmf,~] = cr1nmf( V,K_est ); 
t_inits_cr1nmf = cputime-t0; % cr1-nmf initialization

t0 = cputime; W_spkm = spkm( V,K_est,spkm_maxiter ); 
t_inits_spkm = cputime-t0; % spkm initialization

t0 = cputime;[W_nndsvd,~] = nndsvd( V,K_est );
t_inits_nndsvd = cputime-t0; % nndsvd initialization

% combining with nnlsb
W0 = rand(F,K_est);
[~,~,rel_err_nnlsb,t_nnlsb] = nnlsb(V,W0,nnlsb_maxiter,timelimit); % random initialization
[~,~,rel_err_cr1nmf_nnlsb,t_cr1nmf_nnlsb] = nnlsb(V,W_cr1nmf,nnlsb_maxiter,timelimit-t_inits_cr1nmf); % cr1-nmf initialization
[~,~,rel_err_spkm_nnlsb,t_spkm_nnlsb] = nnlsb(V,W_spkm,nnlsb_maxiter,timelimit-t_inits_spkm); % spkm initialization
[~,~,rel_err_nndsvd_nnlsb,t_nndsvd_nnlsb] = nnlsb(V,W_nndsvd,nnlsb_maxiter,timelimit-t_inits_nndsvd); % nndsvd initialization

nnlsb_mat = [t_nnlsb;rel_err_nnlsb];
cr1nmf_nnlsb_mat = [t_cr1nmf_nnlsb;rel_err_cr1nmf_nnlsb];
spkm_nnlsb_mat = [t_spkm_nnlsb;rel_err_spkm_nnlsb];
nndsvd_nnlsb_mat = [t_nndsvd_nnlsb;rel_err_nndsvd_nnlsb];

t_inits_vec = [t_inits_cr1nmf, t_inits_spkm, t_inits_nndsvd];
filename=['../output/realDataset/t_inits_vec_dataset_',dataset_name,'_K_est',num2str(K_est),'_',str_now,'.mat'];
save(filename,'t_inits_vec');
filename=['../output/realDataset/nnlsb_mat_dataset_',dataset_name,'_K_est',num2str(K_est),'_',str_now,'.mat'];
save(filename,'nnlsb_mat');
filename=['../output/realDataset/cr1nmf_nnlsb_mat_dataset_',dataset_name,'_K_est',num2str(K_est),'_',str_now,'.mat'];
save(filename,'cr1nmf_nnlsb_mat');
filename=['../output/realDataset/spkm_nnlsb_mat_dataset_',dataset_name,'_K_est',num2str(K_est),'_',str_now,'.mat'];
save(filename,'spkm_nnlsb_mat');
filename=['../output/realDataset/nndsvd_nnlsb_mat_dataset_',dataset_name,'_K_est',num2str(K_est),'_',str_now,'.mat'];
save(filename,'nndsvd_nnlsb_mat');

% h = figure('position',[0 0 400 480]);
% plot(nnlsb_mat(1,:), nnlsb_mat(2,:),'r--', 'Linewidth', 1.5); hold on;
% plot(cr1nmf_nnlsb_mat(1,:)+t_inits_vec(1),cr1nmf_nnlsb_mat(2,:),'b', 'Linewidth', 1.5);hold on;
% plot(spkm_nnlsb_mat(1,:)+t_inits_vec(2),spkm_nnlsb_mat(2,:),'g', 'Linewidth', 1.5);hold on;
% plot(nndsvd_nnlsb_mat(1,:)+t_inits_vec(3),nndsvd_nnlsb_mat(2,:),'k-', 'Linewidth', 1.5);hold on;
% xlabel('Running time','FontSize',18);
% ylabel('Relative Error','FontSize',18);
% %title(['Comparison for HALS for ',dataset_name], 'Fontsize',18);
% title(dataset_name, 'Fontsize',18);
% legend('rand','cr1-nmf','spkm','nndsvd','Location', 'Northeast');
% set(gcf,'color','w');str_now=datestr(now,30);set(gca,'FontSize',18);
% filename = ['../Figures/realDataset/relerrTime_NNLSB_dataset_',dataset_name,'_',str_now(1:8),'.pdf'];
% axis([0 300 0.0 0.35])
% export_fig(gcf,'Color','Transparent',filename);
