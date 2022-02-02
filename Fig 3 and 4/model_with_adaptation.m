function [t,xt,x_d,x_r,x_h,g_active] = model_with_adaptation(T,N,K,x0,w_r,w_h,H,bt,g_active_frac,g_type,tau_r,tau_a,sgm,to_plot_V,to_plot_g)
% T - timespan
% N - number of iteration points
% K - number of initial conditions to run
% x0 - initial conditions
% w_r - curvature in x_r
% w_h - curvature in x_h
% H - difference between V(x_r) and V(x_h) (i.e. Delta V)
% bt - coupling strength beta
% g_active_frac - [g_1,g_2] as the locations where the function g is active
% g_type - 1: step function, 2: symmetric triangle, 3: triangle with left angle 90 degree
% tau_r - response of c when growing (rate of increase)
% tau_a - response of c when decreasing 
% sgm - noise strength
% to_plot_V - 1/0 if to plot the potential
% to_plot_g - 1/0 if to plot the function g

%% potential
a = (1/12) * (w_h-w_r) * (w_h+w_r) / H;
b = -(w_h+w_r)^(3/2) / (2*sqrt(3*H));
c = w_r * w_h / (w_h-w_r);

V  = @(x) 1/4*a*x.^4 + 1/3*b*x.^3 + 1/2*c*x.^2;
dV = @(x)     a*x.^3 +     b*x.^2 +     c*x;

x_d = 0;
x_h =  (1/2) * ( -b + sqrt(-4*c*a+b^2) ) / a;
x_r = -(1/2) * (  b + sqrt(-4*c*a+b^2) ) / a;

%% the coupling
h1 = @(x) - x;
h2 = @(x) - x;

g_active = [x_d*(1-g_active_frac(1))+x_r*g_active_frac(1) x_r*(1-g_active_frac(2))+x_h*g_active_frac(2)];  % activity region of function g

switch g_type
    case 1 % simple step
        g = @(x) bt*(x<=g_active(2) & x>=g_active(1));
    case 2 % symmetric triangle
        g = @(x) bt*( max( min(2/(g_active(2)-g_active(1))*(x-g_active(1)) , 2/(g_active(1)-g_active(2))*(x-g_active(2))) , 0 ) ); 
    otherwise % decreasing triangle
        g = @(x) bt * (x<=g_active(2) & x>=g_active(1)) .* ((g_active(2)-x) / (g_active(2)-g_active(1)));
end

%% the model
m = @(t,x) [-dV(x(1,:)) + x(3,:).*h1(x(1,:));
            -dV(x(2,:)) + x(4,:).*h2(x(2,:));
            ( 1/tau_r * (g(x(2,:))-x(3,:)) .* ((g(x(2,:))-x(3,:))>=0) + 1/tau_a * (g(x(2,:))-x(3,:)) .* ((g(x(2,:))-x(3,:))<0) );
            ( 1/tau_r * (g(x(1,:))-x(4,:)) .* ((g(x(1,:))-x(4,:))>=0) + 1/tau_a * (g(x(1,:))-x(4,:)) .* ((g(x(1,:))-x(4,:))<0) )];

sigma_par = [sgm(1);sgm(2);0;0];
if isnan(x0(1))
    x0 = [x_h;x_h;0;0] * ones(1,K);
end

s = @(t,x) reshape([sigma_par(1,:)*(x(1,:)>=0.01);zeros(4,size(x,2));sigma_par(2,:)*(x(2,:)>=0.01);zeros(10,size(x,2))],4,numel(x)); %additive

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
    plot(floor(x_d):.01:ceil(x_h*1.1),g(floor(x_d):.01:ceil(x_h*1.1)),'b','LineWidth',1.3)
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

