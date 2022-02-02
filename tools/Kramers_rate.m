function [r,t] = Kramers_rate(V,c,x_r,x_h,sgm)
%Kramers_rate returns the Kramer's rate and Kramer's time of escape given
%the potential V(x,c), the two equilibria xr (risky/unstable) and xh (helathy/stable)
%and the additive noise strength sgm
r = nan(1,length(c));
t = nan(1,length(c));
for k=1:length(c)
    c_k = c(k);
    x_r_k = x_r(k);
    x_h_k = x_h(k);
    w_r = -MyJacobian(@(x_out) MyJacobian(@(x_in) V(x_in,c_k),x_out,1e-6),x_r_k,1e-3);
    w_h = MyJacobian(@(x_out) MyJacobian(@(x_in) V(x_in,c_k),x_out,1e-6),x_h_k,1e-3);
    H = V(x_r_k,c_k) - V(x_h_k,c_k);
    r(k) = w_r*w_h/2/pi*exp(-2*H/sgm^2);
    t(k) = 1/r(k);
end
end

