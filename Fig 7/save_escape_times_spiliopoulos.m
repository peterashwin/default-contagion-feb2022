close all
clc
clear variables
tic
%% parameters

% randomization
seed = 0;
rng(seed);

% integration parameters
K = 1000; %number of trajectories/initial conditions
T = 50; %timespan
N = T*1e3; %number of timesteps
h = T/N; %time increment
M = 2; %number of dimensions
e = exprnd(1,M,K); %accumulated intensity needed to escape

%parameters of Xt
a = 2;
b = 1;
sigma_0 = 1;
Xt0 = 1;

% parameters of lambda
alph = 4;%*ones(n,1);
lambda_mean = 0.5;%*ones(n,1);
sgm = .9;%*ones(n,1);
Bt = [0 1 2 4]; % <----------------------------------------------------
beta_S = 1;%*ones(n,1);
lambda_0 = .5*ones(M,1);
Xt = cell(K,1);
Wt = cell(K,1);

bt_et_s = struct('bt',cell(1,length(Bt)),'et',cell(1,length(Bt)));

t0=toc;

for bt_iter = 1:length(Bt)
    
    fprintf('\nbeta %g/%g:\n',bt_iter,length(Bt))
    
    %% producing lambda - run until all are defaulted or t=T
    lambda_full = NaN(M,K,N+1);
    Lt_full = NaN(1,K,N+1);
    
    e = exprnd(1,M,K); %accumulated intensity needed to escape
    
    beta_C = Bt(bt_iter);%*ones(n,1);
    for j=1:K
        if ~mod(j,10)
            fprintf('   %g/%g...\n',j,K)
        end
        % dXt = a(b-Xt)dt + sigma*sqrt(Xt)*dWt as in Fig.1 from Spiliopoulos2014
        [~,t,Xt] = MySDE(@(t,x) a*(b-x),@(t,x) sigma_0*sqrt(x),Xt0,[0 T],N);
        
        % generating Wt
        Wt = cumsum([zeros(M,1) sqrt(h)*randn(M,N)],2);
    
        lambda = [lambda_0 NaN(M,N)];
        Lt = [0 NaN(1,N)];
        integral_of_lambda = [zeros(M,1) NaN(M,N)];
        for k=1:N

            if Lt(k) == 1
                break;
            end

            if k == 1
                Lt(k+1) = Lt(k);
            else
                Lt(k+1) = Lt(k) + 1/M * sum(integral_of_lambda(:,k)>e(:,j) & integral_of_lambda(:,k-1)<=e(:,j));
            end
            for l=1:M
                lambda(l,k+1) = lambda(l,k) ...
                    + h * (-alph*(lambda(l,k)-lambda_mean)) ...
                    + sgm*sqrt(lambda(l,k))*(Wt(l,k+1)-Wt(l,k)) ...
                    + beta_C*(Lt(1,k+1)-Lt(1,k)) ...
                    + beta_S*lambda(l,k)*(Xt(1,k+1)-Xt(1,k));
            end
            integral_of_lambda(:,k+1) = integral_of_lambda(:,k) + T/N * lambda(:,k);
        end

        lambda_full(:,j,:) = reshape(lambda,M,1,N+1);
        Lt_full(:,j,:);
        Lt_full(:,j,:) = reshape(Lt,1,1,N+1);
        
        if ~mod(j,100)
            t1 = toc;
            work_done = ((bt_iter-1)+j/K)/length(Bt);
            work_remaining = 1-work_done;
            fprintf('   Estimated time to finish: %ghrs.\n',(t1-t0)/3600*work_remaining/work_done);
        end
        
    end

    %% find escapes on the set of all trajectories
    
    et_s = NaN(M,K);
    for k=1:K
        thr = 1-1/2/M;
        for l=1:M
            et_s(l,k) = EscapeTime(t,1-reshape(Lt_full(1,k,:),1,N+1),thr);
            thr = thr-1/M;
        end
    end

    bt_et_s(bt_iter).bt = beta_C;
    bt_et_s(bt_iter).et = et_s;
end

save([cd() '\saved_escapes\bt_et_s.mat'],'bt_et_s')
% in order to study the escapes in R using copulas, save the output as csv

%%
fprintf('\n');
toc