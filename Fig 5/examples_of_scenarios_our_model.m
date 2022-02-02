close all
clc
clear variables
tic
% run model
T = 20; %time of simulation
N = 1e3*T; %iteration points
K = 1; %number of initial conditions
x0 = NaN; %automatically will put all x in x_h and all c in 0
w_r = 6; w_h = 10; H = 5; %parameters of the potential
bt = 40; % coupling strength 0 or 40
g_active_frac = [.5 .1]; %activation zone of function g
g_type = 2; %1-step, 2-symmetric triangle, 3-triangle with left right angle
tau = 2; %response/adaptation rate
to_plot_V = false; %whether to plot the potential
to_plot_g = to_plot_V; %whether to plot the g function
M = 10; %number of dimensions (M variables x and M variables c)
sgm = 1.5*ones(M,1); %noise strength for x variables (for c it is 0)
if size(sgm,1) > M
    sgm = sgm(1:M);
elseif size(sgm,2) < M
    sgm_corrected = zeros(M,1);
    sgm_corrected(1:size(sgm,1)) = sgm;
    sgm = sgm_corrected;
end
rng(1)

[t,xt,x_d,x_r,x_h,g_active,V] = our_model(T,N,K,M,x0,w_r,w_h,H,bt,g_active_frac,g_type,tau,sgm,to_plot_V,to_plot_g);

%% get one ics and plot it
ic = 1; % 1 <= ic <= K
if K == 1
    xt1 = xt;
else
    xt1 = reshape(xt(:,ic,:),M*2,N+1);
end

lwd1 = {'LineWidth',.7}; lwd2 = {'LineWidth',1.7}; lwd3 = {'LineWidth',1}; fs = {'FontSize',12};

fg1 = figure;
% subplot(10,1,1:4);
hold on;
plot(t(1:50:end),xt1(1:M,1:50:end),lwd1{:});
p(1) =  fill([0 T T 0],g_active([1 1 2 2]),[.9 .9 .9],'LineStyle','none','FaceAlpha',.5,'EdgeAlpha',0);
set(gca,'children',flipud(get(gca,'children')))

plot([0 T],[x_r x_r],'k--')
plot([0 T],[x_h x_h],'k-')

mean_c = mean(xt1(M+1:2*M,:));
x_h_skewed_mean = imag_to_nan(sqrt(3*H)*(w_h+w_r+sqrt((w_h-w_r)^2-4*mean_c*(w_h-w_r)))/(w_h-w_r)/sqrt(w_h+w_r));
x_r_skewed_mean = imag_to_nan(sqrt(3*H)*(w_h+w_r-sqrt((w_h-w_r)^2-4*mean_c*(w_h-w_r)))/(w_h-w_r)/sqrt(w_h+w_r));
x_h_skewed = nan(M,N+1);
x_r_skewed = nan(M,N+1);
for k=1:M
    x_h_skewed(k,:) = imag_to_nan(sqrt(3*H)*(w_h+w_r+sqrt((w_h-w_r)^2-4*xt(M+k,:)*(w_h-w_r)))/(w_h-w_r)/sqrt(w_h+w_r));
    x_r_skewed(k,:) = imag_to_nan(sqrt(3*H)*(w_h+w_r-sqrt((w_h-w_r)^2-4*xt(M+k,:)*(w_h-w_r)))/(w_h-w_r)/sqrt(w_h+w_r));
end
V_skewed = @(x,c) V(x) + 1/2*c.*x.^2;
plot(t,x_h_skewed_mean,'-','Color','r',lwd2{:});
plot(t,x_r_skewed_mean,'--','Color','r',lwd2{:});


text(T+.2,x_h,'x_h',fs{:});
text(T+.2,x_r,'x_r',fs{:});
text(T+.2,x_d,'x_d',fs{:});
box on
xlabel('t','Position',[9 -.2]);ylabel('x','Position',[-.8 2.7],'Rotation',0);
ylim([-.01 6]);
set(gca,'YTick',[0 2 4 6],fs{:});

pos = {'Position',[15 450 700 250]};
ps = {'PaperSize',[16.5 6.7]};
set(fg1,ps{:},pos{:});
print(fg1,[cd() '\figures_examples_of_scenarios\ten_series_x_bt' num2str(bt)],'-dpdf')
% saveas(fg1,[cd() '\figures_examples_of_scenarios\ten_series_x'],'epsc')

%%
if bt>0
    %% c
    fg2 = figure;
    hold on;
    plot(t(1:50:end),xt1(M+1:end,1:50:end),lwd1{:})
    plot(t,mean_c,'r',lwd2{:})
    ylim([0 1.2])
    set(gca,fs{:},'YTick',[0 1.2]);
    box on
    pos = {'Position',[15 60 700 120]};
    ps = {'PaperSize',[16.5 4]};
    xlabel('t','Position',[9 -.03]);ylabel('c','Rotation',0,'Position',[-1 .5]);
    set(fg2,ps{:},pos{:});
    print(fg2,[cd() '\figures_examples_of_scenarios\ten_series_c'],'-dpdf')

    %% H
    fg3 = figure;
    hold on;
    H_skewed_mean = V_skewed(x_r_skewed_mean,mean_c)-V_skewed(x_h_skewed_mean,mean_c);
    H_skewed = V_skewed(x_r_skewed,xt(M+1:end,:)) - V_skewed(x_h_skewed,xt(M+1:end,:));
    plot(t(1:50:end),H_skewed(:,1:50:end),lwd1{:})
    plot(t,H_skewed_mean,'r',lwd2{:})
    xlabel('t','Position',[9 -.2]);ylabel('$\widetilde H$','Rotation',0,'Interpreter','latex','Position',[-1 1.7]);
    set(gca,fs{:});
    xlim([0 T])
    ylim([0 6])
    plot([0 T],[H H],'k--')
    text(T+.2,H,'H',fs{:});
    pos = {'Position',[750 450 700 120]};
    ps = {'PaperSize',[16.5 3.4]};
    set(fg3,ps{:},pos{:});
    box on
    print(fg3,[cd() '\figures_examples_of_scenarios\ten_series_H'],'-dpdf')
end

%% escape rates - ongoging depending on H
if bt>0
    
    subset = 50;
    fg4 = figure;
    hold on;
    kr = NaN(M,size(xt1(:,1:subset:end),2));
    for k=1:M
        fprintf('Finding Kramer''s rate for trajectory %g...\n',k);
        c_k = xt1(M+k,1:subset:end);
        x_h_skewed = imag_to_nan(sqrt(3*H)*(w_h+w_r+sqrt((w_h-w_r)^2-4*c_k*(w_h-w_r)))/(w_h-w_r)/sqrt(w_h+w_r));
        x_r_skewed = imag_to_nan(sqrt(3*H)*(w_h+w_r-sqrt((w_h-w_r)^2-4*c_k*(w_h-w_r)))/(w_h-w_r)/sqrt(w_h+w_r));
        kr(k,:) = Kramers_rate(V_skewed,xt1(M+k,1:subset:end),x_r_skewed,x_h_skewed,sgm(k));
    end

    x_h_skewed = imag_to_nan(sqrt(3*H)*(w_h+w_r+sqrt((w_h-w_r)^2-4*mean_c(1:subset:end)*(w_h-w_r)))/(w_h-w_r)/sqrt(w_h+w_r));
    x_r_skewed = imag_to_nan(sqrt(3*H)*(w_h+w_r-sqrt((w_h-w_r)^2-4*mean_c(1:subset:end)*(w_h-w_r)))/(w_h-w_r)/sqrt(w_h+w_r));
    mean_kr = Kramers_rate(V_skewed,mean_c(1:subset:end),x_r_skewed,x_h_skewed,sgm(k));

    plot(t(1:subset:end),kr(:,:))
    plot(t(1:subset:end),mean_kr,'r',lwd2{:})
    xlabel('t','Position',[9 -.06]);ylabel('r','Rotation',0,'Position',[-.8 .65]);
    ylim([0 1.5]);set(gca,'YTick',[0 1.5]);
    box on
    set(gca,fs{:})
    pos = {'Position',[750 60 700 120]};
    ps = {'PaperSize',[16.5 3.4]};
    set(fg4,ps{:},pos{:});
    print(fg4,[cd() '\figures_examples_of_scenarios\ten_series_r'],'-dpdf')
end
%% write down the escape times
%EscapeTime(t,xt1(1:M,:),x_d/2+x_r/2)

%%
toc