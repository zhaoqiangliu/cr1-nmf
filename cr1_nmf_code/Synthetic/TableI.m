% RUNNING TIMES IN SECONDS OF VARIOUS ALGORITHMS for \alpha = 0.2)
function TableI()
addpath ../algorithms/;
addpath ../plotFuncs/;
clear all;clc;

N_all=10000;
F=1600;K=40;K_est=K;
max_alpha=2*0.2;pert_angle=0.01;

mu_vec=(1:K);
outer_maxit=10;inner_maxit=10;
run_time_mat=zeros(4*inner_maxit*outer_maxit,3);ind_mat=zeros(4*inner_maxit*outer_maxit,3);
timelimit=inf;iter_mult=1000;iter_nnlsb=100;iter_hals=400;
proj_label=1; % take projection
run_time_int=2:4;N_vec=10.^run_time_int; % 10^2, 10^3, 10^4

for outer_iter_num=1:outer_maxit
    fprintf(['The ', num2str(outer_iter_num), 'th outer iteration..\n']);
    [ u_mat,alpha_vec ] = genCones( F,K,max_alpha, pert_angle);
    while(~isempty(find(u_mat(:)<0, 1)))
        fprintf('There are negative values in basis vectors, run again..\n');
        [ u_mat,alpha_vec ] = genCones( F,K,max_alpha, pert_angle);
    end
    
    [ V_all,~ ] = generate_allpts( u_mat,alpha_vec,mu_vec,N_all,proj_label );
    fprintf(['The mean value of length is ',num2str(sum(sum(V_all.*V_all))/N_all),'\n']);
    
    for N=N_vec
        V=V_all(:,1:N);normV=norm(V,'fro');
        for inner_iter_num=1:inner_maxit
            fprintf(['N = ', num2str(N), ', The ', num2str(inner_iter_num), 'th inner iteration..\n']);
            run_time_num=log10(N)-1;t0=cputime;
            [W_cr1nmf,H_cr1nmf] = cr1nmf( V,K_est );
            run_time_mat(inner_iter_num+(outer_iter_num-1)*inner_maxit,run_time_num)=cputime-t0;
            rel_err_cr1nmf=norm(W_cr1nmf*H_cr1nmf-V,'fro')/normV;
            W0=rand(F,K_est);H0=rand(K_est,N);
            [~,~,rel_err_mult,t_mult] = mult(V,W0,H0,iter_mult,timelimit); 
            ind_mult=find(rel_err_mult <= rel_err_cr1nmf, 1 );
            if(~isempty(ind_mult))
                run_time_mat(inner_iter_num+(outer_iter_num-1)*inner_maxit+inner_maxit*outer_maxit,run_time_num)=t_mult(ind_mult);
                ind_mat(inner_iter_num+(outer_iter_num-1)*inner_maxit+inner_maxit*outer_maxit,run_time_num)=ind_mult;
            end
            [~,~,rel_err_nnlsb,t_nnlsb] = nnlsb( V,W0,iter_nnlsb,timelimit ); 
            ind_nnlsb=find(rel_err_nnlsb <= rel_err_cr1nmf, 1 );
            if(~isempty(ind_nnlsb))
                run_time_mat(inner_iter_num+(outer_iter_num-1)*inner_maxit+inner_maxit*outer_maxit*2,run_time_num)=t_nnlsb(ind_nnlsb);
                ind_mat(inner_iter_num+(outer_iter_num-1)*inner_maxit+inner_maxit*outer_maxit*2,run_time_num)=ind_nnlsb;
            end
            [~,~,rel_err_hals,t_hals] = hals(V,W0,H0,iter_hals,timelimit); 
            ind_hals=find(rel_err_hals<=rel_err_cr1nmf, 1 );
            if(~isempty(ind_hals))
                run_time_mat(inner_iter_num+(outer_iter_num-1)*inner_maxit+inner_maxit*outer_maxit*3,run_time_num)=t_hals(ind_hals);
                ind_mat(inner_iter_num+(outer_iter_num-1)*inner_maxit+inner_maxit*outer_maxit*3,run_time_num)=ind_hals;
            end
        end
    end
end

str_now=datestr(now,30);
filename=['../output/Synthetic/run_time_mat_F',num2str(F),'K',num2str(K),'K_est',num2str(K_est),'angle',num2str(50*max_alpha),'_proj_',str_now(1:8),'.mat'];
save(filename,'run_time_mat');
filename=['../output/Synthetic/ind_mat_F',num2str(F),'K',num2str(K),'K_est',num2str(K_est),'angle',num2str(50*max_alpha),'_proj_',str_now(1:8),'.mat'];
save(filename,'ind_mat');
