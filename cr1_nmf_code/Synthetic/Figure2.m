% Errors and performances of various algorithms.
function Figure2(size_angle)
% size_angle = 0.2 or 0.3
addpath ../algorithms/;
addpath ../plotFuncs/;

N_all=10000; % total number of points
F=1600;K=40;K_est=K; % dimensionality and latent dimensionality
max_alpha=2*size_angle;pert_angle=0.01; 
mu_vec=(1:K); % parameters for exponential distributions

alpha_vec = 0.5*max_alpha*ones(1,K);
% non-prob upper bound
nonprob_errbd=max(sin(alpha_vec));
% prob upper bound
f_alpha_vec=sqrt(0.5-sin(2*alpha_vec)./(4*alpha_vec));
prob_errbd=sum(f_alpha_vec.*(mu_vec))/sum(mu_vec);

proj_label=1; % projection to nonnegative orthant
M=100;len=M; 
% we increase N from 10^2 to 10^4 logarithmically
N_vec=sort([round(logspace(2,4,M)),1000]); 
% outer iterations; generate different data matrix (size: F*N_all);
outer_maxit=10;
% inner iterations; 
inner_maxit=10;
% record the results
rel_err_cr1nmf=zeros(inner_maxit*outer_maxit,len+1);rel_err_mult=zeros(inner_maxit*outer_maxit,len+1);
rel_err_nnlsb=zeros(inner_maxit*outer_maxit,len+1);rel_err_hals=zeros(inner_maxit*outer_maxit,len+1);

timelimit=inf;
% number of iterations for mult, nnlsb, hals
iter_mu=100;iter_nnls=20;iter_hals=100;

for outer_iter_num=1:outer_maxit
    fprintf(['The ', num2str(outer_iter_num), 'th outer iteration..\n']);
    [ u_mat,alpha_vec ] = genCones( F,K,max_alpha, pert_angle ); 
    % alpha_vec is fixed as before, but for different outer runs, the u_mat
    % maybe different
    
    % to generate proper parameters u_mat (matrix of basis vectors) and
    % alpha_vec (vector of size angles)
    while(~isempty(find(u_mat(:)<0, 1)))
        fprintf('There are negative values in basis vectors, run again..\n');
        [ u_mat,alpha_vec ] = genCones( F,K,max_alpha, pert_angle );
    end
    
    [ V_all,~ ] = generate_allpts( u_mat,alpha_vec,mu_vec,N_all,proj_label );
    fprintf(['The mean value of length is ',num2str(sum(sum(V_all.*V_all))/N_all),'\n']);
    
    count=0;
    for N = N_vec
        V = V_all(:,1:N);
        normV = norm(V,'fro');
        count = count+1;
        for inner_iter_num = 1:inner_maxit
            fprintf(['N = ', num2str(N), ', The ', num2str(inner_iter_num), 'th inner iteration..\n']);
            [W_cr1nmf,H_cr1nmf] = cr1nmf( V,K_est );
            W0 = rand(F,K_est);H0 = rand(K_est,N);
            [W_mult,H_mult] = mult(V,W0,H0,iter_mu,timelimit);
            [W_nnlsb,H_nnlsb] = nnlsb(V,W0,iter_nnls,timelimit); 
            [W_hals,H_hals] = hals(V,W0,H0,iter_hals,timelimit); 
            rel_err_cr1nmf(inner_iter_num+(outer_iter_num-1)*inner_maxit,count) = norm(W_cr1nmf*H_cr1nmf-V,'fro')/normV;
            rel_err_mult(inner_iter_num+(outer_iter_num-1)*inner_maxit,count) = norm(W_mult*H_mult-V,'fro')/normV;
            rel_err_nnlsb(inner_iter_num+(outer_iter_num-1)*inner_maxit,count) = norm(W_nnlsb*H_nnlsb-V,'fro')/normV;
            rel_err_hals(inner_iter_num+(outer_iter_num-1)*inner_maxit,count) = norm(W_hals*H_hals-V,'fro')/normV;
        end
    end
end

str_now=datestr(now,30);
filename=['../output/Synthetic/rel_err_cr1nmf_F',num2str(F),'K',num2str(K),'K_est',num2str(K_est),'angle',...
    num2str(50*max_alpha),'_proj_',str_now(1:8),'.mat'];
save(filename,'rel_err_cr1nmf');

filename=['../output/Synthetic/rel_err_mult_F',num2str(F),'K',num2str(K),'K_est',num2str(K_est),'angle',...
    num2str(50*max_alpha),'_proj_',str_now(1:8),'.mat'];
save(filename,'rel_err_mult');

filename=['../output/Synthetic/rel_err_nnlsb_F',num2str(F),'K',num2str(K),'K_est',num2str(K_est),'angle',...
    num2str(50*max_alpha),'_proj_',str_now(1:8),'.mat'];
save(filename,'rel_err_nnlsb');

filename=['../output/Synthetic/rel_err_hals_F',num2str(F),'K',num2str(K),'K_est',num2str(K_est),'angle',...
    num2str(50*max_alpha),'_proj_',str_now(1:8),'.mat'];
save(filename,'rel_err_hals');


%% to generate the figures
% % to generate the left figures
% h1 = figure('position',[0 0 450 360]);
% len=length(100:100:N_all);
% semilogx(N_vec, mean(rel_err_cr1nmf), 'm-', 'Linewidth', 1.5); hold on;
% semilogx(N_vec, nonprob_errbd*ones(1,len+1), 'g-.', 'Linewidth', 1.5); hold on;
% semilogx(N_vec, prob_errbd*ones(1,len+1), 'r--', 'Linewidth', 1.5); hold on;
% xlabel('$N$','Interpreter','latex','FontSize',18);
% ylabel('Relative Errors','Interpreter','latex','FontSize',18);
% axis([10^2 10^4 0 0.4])
% [h3,hobj3] = legend('cr1-nmf','non-prob','prob', 'Location', 'Northeast');
% title(['$F=$',num2str(F),', $K=$',num2str(K), ', $\alpha=$',num2str(0.5*max_alpha)],'Interpreter','latex','FontSize',18);
% set(gcf,'color','w');set(gca,'FontSize',18);
% str_now=datestr(now,30);
% filename=['../Figures/Synthetic/checkCorrectness_F',num2str(F),'_K',num2str(K),'K_est',num2str(K_est),'angle',...
%     num2str(50*max_alpha),'_',str_now(1:8),'.pdf'];
% export_fig(gcf,'Color','Transparent',filename);
% 
% % to generate the right figures
% h2 = figure('position',[0 0 450 360]);
% semilogx(N_vec, mean(rel_err_cr1nmf), 'm', 'Linewidth', 1.5); hold on;
% semilogx(N_vec, mean(rel_err_mult), 'b-', 'Linewidth', 1.5, 'MarkerSize',3.0); hold on;
% semilogx(N_vec, mean(rel_err_nnlsb), 'r--', 'Linewidth', 1.5); hold on;
% semilogx(N_vec, mean(rel_err_hals), 'g-.', 'Linewidth', 1.5,'MarkerSize',0.8); hold on;
% xlabel('$N$','Interpreter','latex','FontSize',18);
% %ylabel('Relative Error','Interpreter','latex','FontSize',18);
% title(['$F=$',num2str(F),', $K=$', num2str(K),', $\alpha=$',num2str(0.5*max_alpha)],'interpreter','latex', 'Fontsize',18);
% h4 = legend('cr1-nmf','mult','nnlsb','hals','Location', 'Southeast');
% axis([10^2 10^4 0.08 0.18])
% set(gcf,'color','w');
% set(gca,'FontSize',18);
% str_now=datestr(now,30);
% filename=['../Figures/Synthetic/compareAlgs_F',num2str(F),'K',num2str(K),'K_est',num2str(K_est),'angle',...
%     num2str(50*max_alpha),'_',str_now(1:8),'.pdf'];
% export_fig(gcf,'Color','Transparent',filename);