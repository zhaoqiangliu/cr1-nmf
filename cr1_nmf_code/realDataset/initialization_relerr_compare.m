% this function is for Figure 4 and Table V;

function initialization_relerr_compare(dataset_name)
addpath ../algorithms/;
addpath ../plotFuncs/;

str_now = datestr(now,30);
fprintf(['Now running for dataset ',dataset_name,'...\n']);

% obtaining parameters for the dataset
[ V,K_est,iter_mult,iter_hals,shift_mult,shift_hals,~] = readDataset_relerr_compare(dataset_name);
fprintf('Dataset Parameters are obtained.\n');

max_shift_mult = max(shift_mult); max_shift_hals = max(shift_hals);
[F,N] = size(V); normV = norm(V,'fro');

timelimit = 1e6; % if the running time is too large, terminate the program
maxiter = 10; % the number of runs
nnlsb_maxiter = 100;
% record the results; row number is 4 because we consider 4 types of
% initialization methods, i.e., rand, cr1-nmf, spkm, nndsvd
rel_err_mult_mat = zeros(4,iter_mult-max_shift_mult);
rel_err_hals_mat = zeros(4,iter_hals-max_shift_hals);
rel_err_nnlsb_mat = zeros(4,nnlsb_maxiter); % record the relative error for each iteration for nnlsb algorithm
t_nnlsb_mat = zeros(4*maxiter,nnlsb_maxiter); % record the running time for each iteration for nnlsb algorithm

true_runTime = zeros(maxiter,12); % record results for 4 initialization methods combined with 3 nmf algorithms
t_inits = zeros(maxiter,3); % record running time for 3 initialization methods

rel_err_cr1nmf = zeros(1,maxiter);
spkm_maxiter = 10;
for iter_num = 1:maxiter
    %initializations
    fprintf(['The ',num2str(iter_num),'th iteration..\n']);
    t0 = cputime; [W_cr1nmf,H_cr1nmf] = cr1nmf( V,K_est ); 
    epsilon_cr1nmf = mean(H_cr1nmf(H_cr1nmf>0)); % for perturbation on the right factor matrix
    t_inits(iter_num,1) = cputime-t0; % cr1-nmf initialization
    rel_err_cr1nmf(iter_num) = norm(W_cr1nmf*H_cr1nmf-V,'fro')/normV;
    t0 = cputime; W_spkm = spkm( V,K_est,spkm_maxiter ); H_spkm = rand(K_est,N);  
    t_inits(iter_num,2) = cputime-t0; % spkm initialization
    t0 = cputime; [W_nndsvd,H_nndsvd] = nndsvd( V,K_est ); 
    t_inits(iter_num,3)=cputime-t0; % nndsvd initialization (no need to average over different runs);
    
    %% relative error results
    W0 = rand(F,K_est); H0 = rand(K_est,N);
    % combining with mult
    [~,~,rel_err_mult,t_mult] = mult(V,W0,H0,iter_mult,timelimit); 
    true_runTime(iter_num,1) = t_mult(end);
    rel_err_mult_mat(1,:) = rel_err_mult_mat(1,:)+rel_err_mult((max_shift_mult+1):end); % random initialization
    [~,~,rel_err_cr1nmf_mult,t_cr1nmf_mult] = mult(V,W_cr1nmf,H_cr1nmf+epsilon_cr1nmf*ones(size(H_cr1nmf)),iter_mult-shift_mult(1),timelimit); 
    true_runTime(iter_num,2) = t_cr1nmf_mult(end);
    rel_err_mult_mat(2,:) = rel_err_mult_mat(2,:)+rel_err_cr1nmf_mult((max_shift_mult+1-shift_mult(1)):end); % cr1-nmf initialization
    [~,~,rel_err_spkm_mult,t_spkm_mult] = mult(V,W_spkm,H_spkm,iter_mult-shift_mult(2),timelimit);
    true_runTime(iter_num,3) = t_spkm_mult(end);
    rel_err_mult_mat(3,:) = rel_err_mult_mat(3,:)+rel_err_spkm_mult((max_shift_mult+1-shift_mult(2)):end); % spkm initialization
    [~,~,rel_err_nndsvd_mult,t_nndsvd_mult] = mult(V,W_nndsvd,H_nndsvd,iter_mult-shift_mult(3),timelimit);
    true_runTime(iter_num,4)=t_nndsvd_mult(end);
    rel_err_mult_mat(4,:) = rel_err_mult_mat(4,:)+rel_err_nndsvd_mult((max_shift_mult+1-shift_mult(3)):end); % nndsvd initialization
    
    % combining with nnlsb
    [~,~,rel_err_nnlsb,t_nnlsb] = nnlsb( V,W0,nnlsb_maxiter,timelimit ); 
    true_runTime(iter_num,5)=t_nnlsb(end);
    rel_err_nnlsb_mat(1,:) = rel_err_nnlsb_mat(1,:)+rel_err_nnlsb; t_nnlsb_mat(iter_num,:) = t_nnlsb; % random initialization
    [~,~,rel_err_cr1nmf_nnlsb,t_cr1nmf_nnlsb] = nnlsb( V,W_cr1nmf,nnlsb_maxiter,timelimit );
    true_runTime(iter_num,6)=t_cr1nmf_nnlsb(end);
    rel_err_nnlsb_mat(2,:) = rel_err_nnlsb_mat(2,:)+rel_err_cr1nmf_nnlsb; t_nnlsb_mat(maxiter+iter_num,:) = t_cr1nmf_nnlsb; % cr1-nmf initialization
    [~,~,rel_err_spkm_nnlsb,t_spkm_nnlsb] = nnlsb( V,W_spkm,nnlsb_maxiter,timelimit );
    true_runTime(iter_num,7) = t_spkm_nnlsb(end);
    rel_err_nnlsb_mat(3,:) = rel_err_nnlsb_mat(3,:)+rel_err_spkm_nnlsb; t_nnlsb_mat(maxiter*2+iter_num,:) = t_spkm_nnlsb; % spkm initialization
    [~,~,rel_err_nndsvd_nnlsb,t_nndsvd_nnlsb] = nnlsb( V,W_nndsvd,nnlsb_maxiter,timelimit );
    true_runTime(iter_num,8)=t_nndsvd_nnlsb(end);
    rel_err_nnlsb_mat(4,:) = rel_err_nnlsb_mat(4,:)+rel_err_nndsvd_nnlsb; t_nnlsb_mat(maxiter*3+iter_num,:) = t_nndsvd_nnlsb; % nndsvd initialization
    
    % combining with hals
    [~,~,rel_err_hals,t_hals] = hals(V,W0,H0,iter_hals,timelimit); 
    true_runTime(iter_num,9)=t_hals(end);
    rel_err_hals_mat(1,:) = rel_err_hals_mat(1,:)+rel_err_hals((max_shift_hals+1):end); % random initialization
    [~,~,rel_err_cr1nmf_hals,t_cr1nmf_hals] = hals(V,W_cr1nmf,H_cr1nmf+epsilon_cr1nmf*ones(size(H_cr1nmf)),iter_hals-shift_hals(1),timelimit);
    true_runTime(iter_num,10)=t_cr1nmf_hals(end);
    rel_err_hals_mat(2,:) = rel_err_hals_mat(2,:)+rel_err_cr1nmf_hals((max_shift_hals+1-shift_hals(1)):end); % cr1-nmf initialization
    [~,~,rel_err_spkm_hals,t_spkm_hals] = hals(V,W_spkm,H_spkm,iter_hals-shift_hals(2),timelimit);
    true_runTime(iter_num,11)=t_spkm_hals(end);
    rel_err_hals_mat(3,:) = rel_err_hals_mat(3,:)+rel_err_spkm_hals((max_shift_hals+1-shift_hals(2)):end); % spkm initialization
    [~,~,rel_err_nndsvd_hals,t_nndsvd_hals] = hals(V,W_nndsvd,H_nndsvd,iter_hals-shift_hals(3),timelimit);
    true_runTime(iter_num,12)=t_nndsvd_hals(end);
    rel_err_hals_mat(4,:) = rel_err_hals_mat(4,:)+rel_err_nndsvd_hals((max_shift_hals+1-shift_hals(3)):end); % nndsvd initialization
end
rel_err_mult_mat = rel_err_mult_mat/maxiter; 
rel_err_nnlsb_mat = rel_err_nnlsb_mat/maxiter; 
rel_err_hals_mat = rel_err_hals_mat/maxiter;

filename = ['../output/realDataset/rel_err_mult_mat_dataset_',dataset_name,'_K_est',num2str(K_est),'_',str_now(1:8),'.mat'];
save(filename,'rel_err_mult_mat');
filename = ['../output/realDataset/rel_err_nnlsb_mat_dataset_',dataset_name,'_K_est',num2str(K_est),'_',str_now(1:8),'.mat'];
save(filename,'rel_err_nnlsb_mat');
filename = ['../output/realDataset/rel_err_hals_mat_dataset_',dataset_name,'_K_est',num2str(K_est),'_',str_now(1:8),'.mat'];
save(filename,'rel_err_hals_mat');
filename = ['../output/realDataset/rel_err_cr1nmf_dataset_',dataset_name,'_K_est',num2str(K_est),'_',str_now(1:8),'.mat'];
save(filename,'rel_err_cr1nmf');
filename = ['../output/realDataset/true_runTime_dataset_',dataset_name,'_K_est',num2str(K_est),'_',str_now(1:8),'.mat'];
save(filename,'true_runTime');
filename = ['../output/realDataset/t_inits_dataset_',dataset_name,'_K_est',num2str(K_est),'_',str_now(1:8),'.mat'];
save(filename,'t_inits');
filename = ['../output/realDataset/t_nnlsb_mat_dataset_',dataset_name,'_K_est',num2str(K_est),'_',str_now(1:8),'.mat'];
save(filename,'t_nnlsb_mat');


%%
% str_now=datestr(now,30);normV=norm(V,'fro');
% min_relerrmu=min(rel_err_mult_mat(:,end))/normV;max_relerrmu=max(rel_err_mult_mat(:,1))/normV;
% errmu_low=max([0.5*min_relerrmu,min_relerrmu-0.01]);
% errmu_high=min([max_relerrmu+0.01,2*max_relerrmu]);
% h1 = figure('position',[0 0 300 360]);
% plot(rel_err_mult_mat(1,:),'r--', 'Linewidth', 1.5);hold on;
% plot(rel_err_mult_mat(2,:),'b', 'Linewidth', 1.5);hold on;
% plot(rel_err_mult_mat(3,:),'g', 'Linewidth', 1.5);hold on;
% plot(rel_err_mult_mat(4,:),'k-.', 'Linewidth', 1.5);hold on;
% xlabel('Iteration Number','FontSize',18);
% ylabel('Relative Error','FontSize',18);
% title('MULT', 'Fontsize',18);
% legend('rand','cr1-nmf','spkm','nndsvd','Location', 'Northeast');
% set(gcf,'color','w');set(gca,'FontSize',18);
% filename=['../Figures/realDataset/MULT_dataset_',dataset_name,'_K_est',num2str(K_est),'_',str_now(1:8),'.pdf'];
% axis([1 size(rel_err_mult_mat,2) errmu_low errmu_high])
% export_fig(gcf,'Color','Transparent',filename);
% 
% 
% min_relerrnnls = min(rel_err_nnlsb_mat(:,end))/normV;max_relerrnnls=max(rel_err_nnlsb_mat(:,1))/normV;
% errnnls_low = max([0.5*min_relerrnnls,min_relerrnnls-0.01]);
% errnnls_high = min([max_relerrnnls+0.01,2*max_relerrnnls]);
% h2 = figure('position',[0 0 300 360]);
% plot(rel_err_nnlsb_mat(1,:),'r--', 'Linewidth', 1.5);hold on;
% plot(rel_err_nnlsb_mat(2,:),'b', 'Linewidth', 1.5);hold on;
% plot(rel_err_nnlsb_mat(3,:),'g', 'Linewidth', 1.5);hold on;
% plot(rel_err_nnlsb_mat(4,:),'k-', 'Linewidth', 1.5);hold on;
% %xlabel('Iteration Number','FontSize',18);
% %ylabel('Relative Error','FontSize',18);
% %title('NNLSB', 'Fontsize',18);
% legend('rand','cr1-nmf','spkm','nndsvd','Location', 'Northeast');
% set(gcf,'color','w');str_now=datestr(now,30);set(gca,'FontSize',18);
% filename = ['../Figures/realDataset/NNLSB_dataset_',dataset_name,'_K_est',num2str(K_est),'_',str_now(1:8),'.pdf'];
% axis([1 size(rel_err_nnlsb_mat,2) errnnls_low errnnls_high])
% export_fig(gcf,'Color','Transparent',filename);
% 
% 
% min_relerrhals = min(rel_err_hals_mat(:,end))/normV;max_relerrhals=max(rel_err_hals_mat(:,1))/normV;
% errhals_low = max([0.5*min_relerrhals,min_relerrhals-0.01]);
% errhals_high = min([max_relerrhals+0.01,2*max_relerrhals]);
% h3 = figure('position',[0 0 300 360]);
% plot(rel_err_hals_mat(1,:),'r--', 'Linewidth', 1.5);hold on;
% plot(rel_err_hals_mat(2,:),'b', 'Linewidth', 1.5);hold on;
% plot(rel_err_hals_mat(3,:),'g', 'Linewidth', 1.5);hold on;
% plot(rel_err_hals_mat(4,:),'k-', 'Linewidth', 1.5);hold on;
% %xlabel('Iteration Number','FontSize',18);
% %ylabel('Relative Error','FontSize',18);
% %title('HALS', 'Fontsize',18);
% legend('rand','cr1-nmf','spkm','nndsvd','Location', 'Northeast');
% set(gcf,'color','w');str_now=datestr(now,30);set(gca,'FontSize',18);
% filename = ['../Figures/realDataset/HALS_dataset_',dataset_name,'_K_est',num2str(K_est),'_',str_now(1:8),'.pdf'];
% axis([1 size(rel_err_hals_mat,2) errhals_low errhals_high])
% export_fig(gcf,'Color','Transparent',filename);
% 
% 
