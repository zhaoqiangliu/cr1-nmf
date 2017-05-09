% Estimated number of circular cones K with different noise levels.
function Figure3(noise_level)
% noise level = 0.1 or 0.5
addpath ../algorithms/;
addpath ../plotFuncs/;

%clear all;clc;
size_angle=0.3;
N_all=10000;
F=1600;K=40;
max_alpha=2*size_angle;pert_angle=0.01;
fprintf(['The size angle = ', num2str(size_angle), '..\n']);

mu_vec=(1:K);
proj_label=1;% take projection
M=100;len=M;
N_vec=sort([round(logspace(2,4,M)),1000]);
outer_maxit=1000;
estK_mat=zeros(outer_maxit,len+1);
%noise_level=0.1;
K_min=11;K_max=200;

t0=cputime;
for outer_iter_num=1:outer_maxit
    fprintf(['The ', num2str(outer_iter_num), 'th outer iteration..\n']);
    [ u_mat,alpha_vec ] = genCones( F,K,max_alpha, pert_angle);
    while(~isempty(find(u_mat(:)<0, 1)))
        fprintf('There are negative values in basis vectors, run again..\n');
        [ u_mat,alpha_vec ] = genCones( F,K,max_alpha, pert_angle);
    end
    
    [ V_all,~ ] = generate_allpts( u_mat,alpha_vec,mu_vec,N_all,proj_label );
    fprintf(['The mean value of length is ',num2str(sum(sum(V_all.*V_all))/N_all),'\n']);
    V_all=V_all+noise_level*randn(F,N_all); % perturbation
    length(find(max(V_all)<0))
    V_all=V_all.*(V_all>0);
    
    count=0;
    for N=N_vec
        fprintf(['The ', num2str(outer_iter_num), 'th outer iteration...Now we consider N = ',num2str(N),'...\n']);
        V=V_all(:,1:N);count=count+1;
        [~,Sigma]=mySVD(V,min(K_max+1,N));
        diagS=diag(Sigma)';diagS=diagS(K_min:end);
        sv_ratio=diagS(1:(end-1))./diagS(2:end);
        [~,ind]=max(sv_ratio);
        estK_mat(outer_iter_num,count)=ind+K_min-1;
    end
end
total_time=cputime-t0;

% save to file
str_now=datestr(now,30);
filename=['../output/Synthetic/estK_mat_F',num2str(F),'K',num2str(K),'_noiselevel',num2str(noise_level*10),'_angle',...
    num2str(50*max_alpha),'_proj_',str_now(1:8),'.mat'];
save(filename,'estK_mat');


% % plot the results
% M=100;
% h1 = figure('position',[0 0 450 360]);
% len=length(100:100:N_all);
% ha=errorbar(N_vec,mean(estK_mat),std(estK_mat),'b-','MarkerSize',1.6,'MarkerEdgeColor','m','MarkerFaceColor','m'); hold on;
% %errorbar_tick(ha,20);
% hb = get(ha,'children'); 
% Xdata = get(hb(2),'Xdata');
% temp = 4:3:length(Xdata);
% temp(3:3:end) = [];
% % xleft and xright contain the indices of the left and right endpoints of the horizontal lines
% xleft = temp; xright = temp+1;
% % Increase line length by 0.2 units
% Xdata(xleft) = log(Xdata(xleft) +10);
% Xdata(xright) = log(Xdata(xright) -10);
% plot(N_vec, K*ones(1,len+1), 'r--', 'Linewidth', 1.5); hold on;
% xlabel('$N$','Interpreter','latex','FontSize',20);
% ylabel('$K$','Interpreter','latex','FontSize',20);
% set(gca,'xscale','log')
% set(hb(2),'Xdata',Xdata)
% axis([10^2 10^4 0 K_max])
% [h,hobj] = legend('$\hat{K}$','true $K$', 'Location', 'Northeast');
% set(h, 'Interpreter', 'latex')
% title(['$F=$ ',num2str(F),', true $K=$ ',num2str(K), ', $\alpha=$ ',num2str(0.5*max_alpha),', $\delta=$ ',num2str(noise_level)],...
%     'Interpreter','latex','FontSize',18);
% set(gcf,'color','w');set(gca,'FontSize',18);
% str_now=datestr(now,30);
% filename=['../Figures/Synthetic/deterK_F',num2str(F),'_K',num2str(K),'_noiselevel',num2str(noise_level*10),'_angle',...
%     num2str(50*max_alpha),'_',str_now(1:8),'.pdf'];
% export_fig(gcf,'Color','Transparent',filename);
% 
