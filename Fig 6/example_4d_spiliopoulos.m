close all
clc
clear variables
tic
datetime('now')
%% parameters

% randomization
seed = 12;
rng(seed);

% integration parameters
K = 1; %number of trajectories/initial conditions
T = 6;%35; %timespan
N = T*1e3; %number of timesteps
h = T/N; %time increment
n = 4; %number of dimensions
e = exprnd(1,n,K); %accumulated intensity needed to escape
%% generating Xt
% dXt = a(b-Xt)dt + sigma*sqrt(Xt)*dWt as in Fig.1 from Spiliopoulos2014
a = 2;
b = 1;
sigma_0 = 1;
Xt0 = ones(1,K);
[~,t,Xt] = MySDE(@(t,x) a*(b-x),@(t,x) sigma_0*sqrt(x),Xt0,[0 T],N);

%% generating Wt
Wt = cumsum(cat(3,zeros(n,K,1),sqrt(h)*randn(n,K,N)),3);

%% defining Lt
% L = @(tt,e,lambda) 1/n*sum((sum(lambda(:,t<=tt),2)>=e));
%tt - current time in simulation
%e - exponential r.v.
%% producing lambda
lambda_full = NaN(n,N+1);
Lt_full = NaN(1,N+1);

alph = 4;%*ones(n,1);
lambda_mean = 0.5;%*ones(n,1);
sgm = .9;%*ones(n,1);
beta_C_val = 2;
beta_C = beta_C_val;%*ones(n,1);
beta_S = 1;%*ones(n,1);
lambda_0 = .5*ones(n,1);


lambda = [lambda_0 NaN(n,N)];
Lt = [0 NaN(1,N)];
integral_of_lambda = [zeros(n,1) NaN(n,N)];
for k=1:N
    if k == 1
        Lt(k+1) = Lt(k);
    else
        Lt(k+1) = Lt(k) + 1/n * sum(integral_of_lambda(:,k)>e(:) & integral_of_lambda(:,k-1)<=e(:));
    end
    for l=1:n
        lambda(l,k+1) = lambda(l,k) ...
            + h * (-alph*(lambda(l,k)-lambda_mean)) ...
            + sgm*sqrt(lambda(l,k))*(Wt(l,k+1)-Wt(l,k)) ...
            + beta_C*(Lt(1,k+1)-Lt(1,k)) ...
            + beta_S*lambda(l,k)*(Xt(1,k+1)-Xt(1,k));
    end
    integral_of_lambda(:,k+1) = integral_of_lambda(:,k) + T/N * lambda(:,k);
end

lambda_full(:,:) = reshape(lambda,n,1,N+1);
Lt_full(:,:);
Lt_full(:,:) = reshape(Lt,1,1,N+1);

lambda_to_plot = lambda;
integral_of_lambda_to_plot = integral_of_lambda;
for k=1:n
    lambda_to_plot(k,integral_of_lambda(k,:)>e(k)) = nan;
    integral_of_lambda_to_plot(k,integral_of_lambda(k,:)>e(k)) = nan;
end

lwd = {'LineWidth',1};
fg=figure('Position',[50,50,750,350],'PaperSize',[19.5 9.5]);hold on;
set(gca,'FontSize',15)
col = 'rgbc';
for k=1:n
    plot(t,lambda_to_plot(k,:),'Color',col(k),lwd{:})
end
ylim([0 2.5])
xlabel('time');ylabel('intensity')
yyaxis right
plot(t,Lt,'k--',lwd{:})
ylim([0 1.1])
ylabel('loss')
set(gca,'YColor','k')
box on
n_round = 2;
legend(['\lambda_1, \epsilon_1=' num2str(round(e(1),n_round))],...
       ['\lambda_2, \epsilon_2=' num2str(round(e(2),n_round)) '0'],... %added zero here!
       ['\lambda_3, \epsilon_3=' num2str(round(e(3),n_round))],...
       ['\lambda_4, \epsilon_4=' num2str(round(e(4),n_round))],...
       'L_t','Location','East')
print(fg,[cd() '\figures\ex_spil'],'-dpdf')

%%
toc