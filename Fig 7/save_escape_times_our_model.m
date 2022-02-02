close all
clc
clear variables
% run model
T = 20; %time of simulation
N = 1e3*T; %iteration points
K = 1000; %number of initial conditions
x0 = NaN; %automatically will put all x in x_h and all c in 0
w_r = 6; w_h = 10; H = 5; %parameters of the potential
g_active_frac = [.5 .1]; %activation zone of function g
g_type = 2; %1-step, 2-symmetric triangle, 3-triangle with left right angle
tau = 2; %response/adaptation rate
to_plot_V = false; %whether to plot the potential
to_plot_g = to_plot_V; %whether to plot the g function

M = 2; %number of dimensions (M variables x and M variables c)

sgm = 1.5*ones(M,1); %noise strength for x variables (for c it is 0)
if size(sgm,1) > M
    sgm = sgm(1:M);
elseif size(sgm,2) < M
    sgm_corrected = zeros(M,1);
    sgm_corrected(1:size(sgm,1)) = sgm;
    sgm = sgm_corrected;
end
rng(0)

%%
tic;t0 = toc;
Bt = [0 2.5 5 10]; % <----------------------------------------------------
all_data = cell(length(Bt),2);

Bt_iter = [1 2 3 4];

for bt_iter=Bt_iter
    
    bt_et_partitioned = struct('tau',cell(1,1),'bt',cell(1,1),'et',cell(1,1));
    bt = Bt(bt_iter);
    fprintf('\nbeta %g/%g\n',bt_iter,length(Bt))
    [t,xct,x_d,x_r,x_h,g_active,V] = our_model(T,N,K,M,x0,w_r,w_h,H,bt,g_active_frac,g_type,tau,sgm,to_plot_V,to_plot_g);
    xt = xct(1:M,:,:);
    ct = xct(M+1:end,:,:);

    % find the escape times
    et = EscapeTime(t,xt,x_d/2+x_r/2);
    all_traj = struct('t',cell(1,K),'x',cell(1,K),'c',cell(1,K),'et',cell(1,K));

    % if there is nan in any of trajectories, continue for time T
    for k=1:K
        fprintf('\nbeta %g/%g, ics %g/%g\n',bt_iter,length(Bt),k,K)
        xtk = reshape(xt(:,k,:),M,N+1);
        ctk = reshape(ct(:,k,:),M,N+1);
        tk = t;
        l = 0;
        while any(isnan(et(:,k)))
            l = l+1;
            fprintf('beta %g/%g, ics %g/%g integrating further (%g)\n',bt_iter,length(Bt),k,K,l)
            [t_new,xct_new] = our_model(T,N,1,M,[xtk(:,end);ctk(:,end)],w_r,w_h,H,bt,g_active_frac,g_type,tau,sgm,to_plot_V,to_plot_g);
            tk = [tk(1:end-1) t_new+l*T];
            xtk = [xtk(:,1:end-1) xct_new(1:M,:)];
            ctk = [ctk(:,1:end-1) xct_new(M+1:end,:)];
            et(:,k) = EscapeTime(tk,xtk,x_d/2+x_r/2);
        end
        all_traj(k).t = tk;
        all_traj(k).x = xtk;
        all_traj(k).c = ctk;
        all_traj(k).et = et(:,k);
        
        t1 = toc;
        work_done = ((find(bt_iter==Bt_iter)-1)+k/K)/length(Bt_iter);
        work_remaining = 1-work_done;
        fprintf('Estimated time to finish: %ghrs.\n',(t1-t0)/3600*work_remaining/work_done);
    end
    
    bt_et_partitioned.tau = tau;
    bt_et_partitioned.bt = bt;
    bt_et_partitioned.et = sort(et,1);
    
    save([cd() '\saved_escapes\partitioned\bt_et' num2str(bt_iter) '.mat'],'bt_et_partitioned')
    
%     clear bt_et_partitioned
    t1 = toc;
end


%% load the partitions and save as one file
bt_et = struct('tau',cell(1,length(Bt)),'bt',cell(1,length(Bt)),'et',cell(1,length(Bt)));
for bt_iter=1:length(Bt)
    load([cd() '\saved_escapes\partitioned\bt_et' num2str(bt_iter) '.mat']);
    bt_et(bt_iter) = bt_et_partitioned;
end
save([cd() '\saved_escapes\bt_et.mat'],'bt_et')
% in order to study the escapes in R using copulas, save the output as csv

%%

run save_escape_times_5d_spiliopoulos.m


toc
