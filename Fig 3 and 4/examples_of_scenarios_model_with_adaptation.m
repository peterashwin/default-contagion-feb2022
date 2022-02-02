% close all
clc
clear variables
tic
% run model
T = 10;
N = 1e4;
K = 1;
x0 = NaN;
w_r = 6; w_h = 10; H = 5;
bt = 16; %8,12,16
g_active_frac = [.5 .1];
g_type = 2;
tau_r = 2;
tau_a = tau_r;
sgm = [2 0];
to_plot_V = true;
to_plot_g = true;
rng(1)

[t,xt,x_d,x_r,x_h,g_active] = model_with_adaptation(T,N,K,x0,w_r,w_h,H,bt,g_active_frac,g_type,tau_r,tau_a,sgm,to_plot_V,to_plot_g);
fg = figure(1);
annotation(fg,'textarrow',[0.397 0.327],[0.289 0.152],'String',{'x_{dr}'},'FontSize',12);
annotation(fg,'textarrow',[0.60 0.538],[0.269 0.149],'String',{'x_{rh}'},'FontSize',12)
set(fg,'PaperSize',[14 5.5]);
warning('off')
print(fg,[cd() '\figures_examples_of_scenarios\gfunction'],'-dpdf')
warning('on')

min_stable   = sqrt(3*H)*(w_h+w_r+sqrt((w_h-w_r)^2-4*bt*(w_h-w_r)))/(w_h-w_r)/sqrt(w_h+w_r);
max_unstable = sqrt(3*H)*(w_h+w_r-sqrt((w_h-w_r)^2-4*bt*(w_h-w_r)))/(w_h-w_r)/sqrt(w_h+w_r);
disp(min_stable);disp(max_unstable);

%% get one trajectory and plot it
ic = 1; % 1 <= ic <= K
if K == 1
    xt1 = xt;
else
    xt1 = reshape(xt(:,ic,:),4,N+1);
end

lwd = {'LineWidth',1.3}; fs = {'FontSize',14};

fg = figure;
% subplot(2,1,1);
hold on;
p(3) =  fill([0 T T 0],g_active([1 1 2 2]),[.9 .9 .9],'LineStyle','none','FaceAlpha',.5,'EdgeAlpha',0);
plot([0 T],[x_r x_r],'k--')
plot([0 T],[x_h x_h],'k-')
% set(gca,'children',flipud(get(gca,'children')),fs{:})
% axis([0 T 0 ceil(max(xt1(:)))])

%how far can equilibrium go:
min_stable   = sqrt(3*H)*(w_h+w_r+sqrt((w_h-w_r)^2-4*bt*(w_h-w_r)))/(w_h-w_r)/sqrt(w_h+w_r);
max_unstable = sqrt(3*H)*(w_h+w_r-sqrt((w_h-w_r)^2-4*bt*(w_h-w_r)))/(w_h-w_r)/sqrt(w_h+w_r);
disp(min_stable);disp(max_unstable);

%shifted red equilibria
x_h_skewed = imag_to_nan(sqrt(3*H)*(w_h+w_r+sqrt((w_h-w_r)^2-4*xt1(4,:)*(w_h-w_r)))/(w_h-w_r)/sqrt(w_h+w_r));
x_r_skewed = imag_to_nan(sqrt(3*H)*(w_h+w_r-sqrt((w_h-w_r)^2-4*xt1(4,:)*(w_h-w_r)))/(w_h-w_r)/sqrt(w_h+w_r));
p(4) = plot(t,x_h_skewed,'-','Color',[1 0 0]);
p(5) = plot(t,x_r_skewed,'--','Color',[1 0 0]);
%shifted blue equilibria
x_h_skewed = imag_to_nan(sqrt(3*H)*(w_h+w_r+sqrt((w_h-w_r)^2-4*xt1(3,:)*(w_h-w_r)))/(w_h-w_r)/sqrt(w_h+w_r));
x_r_skewed = imag_to_nan(sqrt(3*H)*(w_h+w_r-sqrt((w_h-w_r)^2-4*xt1(3,:)*(w_h-w_r)))/(w_h-w_r)/sqrt(w_h+w_r));
p(6) = plot(t,x_h_skewed,'-','Color',[0 0 1]);
p(7) = plot(t,x_r_skewed,'--','Color',[0 0 1]);

p(1) = plot(t,xt1(1,:),'b',lwd{:});
p(2) = plot(t,xt1(2,:),'r',lwd{:});

text(x_h,T,'x_h');
text(x_r,T,'x_r');
text(x_d,T,'x_d');
box on
xlabel('t','Position',[7.4 -.2]);ylabel('x','Position',[-.4 4.9],'Rotation',0);
ylim([-.01 6]);
set(gca,fs{:},'YTick',[0 2 4 6]);
% h_lgd = legend(p(:),'x_1','x_2','activation of g function','real stable eq. for x_1','real unstable eq. for x_1',...
%        'real stable eq. for x_2','real unstable eq. for x_2');
% set(h_lgd,'FontSize',8,'Location','SouthWest');
if bt == 8
    md = '1';
elseif bt == 12
    md = '2';
elseif bt == 16
    md = '3';
end

pos = {'Position',[800 250 350 250]};
ps = {'PaperSize',[8 7]};
set(fg,ps{:},pos{:});
print(fg,[cd() '\figures_examples_of_scenarios\x_mode' md],'-dpdf')

% subplot(2,1,2);
fg = figure;
hold on;
plot(t,xt1(3,:),'b',lwd{:})
plot(t,xt1(4,:),'r',lwd{:})
set(gca,fs{:},'YTick',[0 1 2])
set(fg,ps{:},pos{:})
xlabel('t','Position',[7.4 -.065]);ylabel('c','Position',[-.4 1.5],'Rotation',0);
ylim([0 2])
box on
print(fg,[cd() '\figures_examples_of_scenarios\c_mode' md],'-dpdf')

%%
toc