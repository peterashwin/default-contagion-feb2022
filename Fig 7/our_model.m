function [t,xt,x_d,x_r,x_h,g_active,V] = our_model(T,N,K,M,x0,w_r,w_h,H,bt,g_active_frac,g_type,tau,sgm,to_plot_V,to_plot_g)
% T - timespan
% N - number of iteration points
% K - number of initial conditions to run
% M - how many players to include
% x0 - initial conditions
% w_r - curvature in x_r
% w_h - curvature in x_h
% H - difference between V(x_r) and V(x_h) (i.e. Delta V)
% bt - coupling strength beta
% g_active_frac - [g_1,g_2] as the locations where the function g is active
% g_type - 1: step function, 2: symmetric triangle, 3: triangle with left angle 90 degree
% tau - response of c (=tau_r=tau_a)
% sgm - noise strength
% to_plot_V - 1/0 if to plot the potential
% to_plot_g - 1/0 if to plot the function g

%% potential
a = (1/12) * (w_h-w_r) * (w_h+w_r) / H;
b = -(w_h+w_r)^(3/2) / (2*sqrt(3*H));
c = w_r * w_h / (w_h-w_r);

V   = @(x) 1/4*a*x.^4 + 1/3*b*x.^3 + 1/2*c*x.^2;
dV  = @(x)     a*x.^3 +     b*x.^2 +     c*x;
% ddV = @(x)   3*a*x^2  +   2*b*x    +     c;

x_d = 0;
x_h =  (1/2) * ( -b + sqrt(-4*c*a+b^2) ) / a;
x_r = -(1/2) * (  b + sqrt(-4*c*a+b^2) ) / a;

%% the coupling
h = @(x) - x;

g_active = [x_d*(1-g_active_frac(1))+x_r*g_active_frac(1) x_r*(1-g_active_frac(2))+x_h*g_active_frac(2)];  % activity region of function g

switch g_type
    case 1 % simple step
        g = @(x) bt / (M-1) * (x<=g_active(2) & x>=g_active(1));
    case 2 % symmetric triangle
        g = @(x) bt / (M-1) * ( max( min(2/(g_active(2)-g_active(1))*(x-g_active(1)) , 2/(g_active(1)-g_active(2))*(x-g_active(2))) , 0 ) ); 
    otherwise % decreasing triangle
        g = @(x) bt / (M-1) * (x<=g_active(2) & x>=g_active(1)) .* ((g_active(2)-x) / (g_active(2)-g_active(1)));
end

%% the model

% subfunction_x
fileID = fopen('subfunction_x.m','w');
fprintf(fileID,'function y = subfunction_x(x,dV,h)\n');
fprintf(fileID,'y = [...\n');
for M_cnt = 1:M
    fprintf(fileID,['-dV(x(' num2str(M_cnt) ',:)) + x(' num2str(M_cnt+M) ',:).*h(x(' num2str(M_cnt) ',:));...\n']);
end
fprintf(fileID,'];\n');
fprintf(fileID,'end');
fclose(fileID);

% subfunction_c
fileID = fopen('subfunction_c.m','w');
fprintf(fileID,'function y = subfunction_c(x,tau,g)\n');
fprintf(fileID,'y = [...\n');
for M_cnt1 = 1:M
    sum_g_txt = '';
    for M_cnt2 = 1:M
        if M_cnt2 ~= M_cnt1
            sum_g_txt = [sum_g_txt 'g(x(' num2str(M_cnt2) ',:))+'];
        end
    end
    sum_g_txt = sum_g_txt(1:end-1);
    fprintf(fileID,['1/tau * (' sum_g_txt '-x(' num2str(M_cnt1+M) ',:));...\n']);
end
fprintf(fileID,'];\n');
fprintf(fileID,'end');
fclose(fileID);

m = @(t,x) [subfunction_x(x , dV , h);
            subfunction_c(x , tau , g)];

sigma_par = [sgm;zeros(size(sgm))];
if isnan(x0(1))
    x0 = [x_h*ones(M,1);zeros(M,1)] * ones(1,K);
end

s = @(t,x) diag_mat(repmat(sigma_par,1,size(x,2)).*(x>=0.01)); %additive

[~,t,xt] = MySDE(m,s,x0,[0 T],N);

%% plots
if to_plot_V && ~to_plot_g % only V
    figure('Position',[50 500 300 100])
    hold on
    plot(floor(x_d):.01:ceil(x_h*1.1),V(floor(x_d):.01:ceil(x_h*1.1)),'k','LineWidth',1.3)
    xlabel('x');ylabel('V(x)')
    box on
    set(gca,'FontSize',12)
elseif ~to_plot_V && to_plot_g % only g
    figure('Position',[50 500 300 100])
    hold on
    plot(floor(x_d):.01:ceil(x_h*1.1),g(floor(x_d):.01:ceil(x_h*1.1)),'b','LineWidth',1.3)
    ylim([0 1.05])
    xlabel('x');ylabel('g(x)')
    box on
    set(gca,'FontSize',12)
elseif to_plot_V && to_plot_g % both V and g
    figure('Position',[50 500 600 200])
    hold on
    yyaxis left;set(gca,'YColor','k')
    plot(floor(x_d):.01:ceil(x_h*1.1),V(floor(x_d):.01:ceil(x_h*1.1)),'k','LineWidth',1.3)
    text(x_r+0.05,V(x_r)-1.5,'x_r','FontSize',12);
    text(x_h+0.05,V(x_h)-1.5,'x_h','FontSize',12);
    text(x_d+0.05,V(x_d)+2,'x_d','FontSize',12);
    yyaxis right;set(gca,'YColor','b')
    plot(floor(x_d):.01:ceil(x_h*1.1),g(floor(x_d):.01:ceil(x_h*1.1))*(M-1),'b','LineWidth',1.3)
    plot([x_r x_r],[0 bt+.05],'k--');plot([x_h x_h],[0 bt+.05],'k--');
    ylim([0 bt+.05])
    legend('V(x)','g(x)','Location','NorthWest')
    box on
    set(gca,'FontSize',12)
    ytck_right = get(gca,'YTick');
    set(gca,'YTick',ytck_right([1 end]));
    yt_label = get(gca,'YTickLabel');
    yt_label{end} = '1';
    set(gca,'YTickLabel',yt_label);
%     xlabel('x');yyaxis left;ylabel('V(x)');yyaxis right;ylabel('g(x)');
    if ~mod(ceil(x_h*1.1),2)
        xtck = get(gca,'XTick');
        xtck = xtck(1:2:end);
        set(gca,'XTick',xtck);
    end
    xlbl = xlabel('x');
    x_pos = get(xlbl,'Position') + [0 3 0];
    set(xlbl,'Position',x_pos);
end

end

