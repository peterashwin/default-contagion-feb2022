function [xend,t,xt] = MySDE(m,s,x0,tspan,N,varargin)
%% [xend,t,xt] = MySDE(m,s,x0,tspan,N)
%   Inputs:
%       m: function defining the deterministic part of SDE.
%          m should accept two arguments: t (a scalar) and
%          x (an n-dimensional vector). The function m should
%          return an n-dimensional vector y (the time derivative).
%          Typical calling sequence: y=m(t,x), returning the value
%          of m at time t in position x.
%          m can be also defined to work on multiple values of x - it is
%          much quicker then for multiple initial conditions and in that
%          case m returns a matrix of the same size as x0.
%       s: function defining the stochastic part of SDE. s should
%          accept two arguments: t (a scalar) and x (an
%          n-dimensional vector). The function s should return an
%          n x n-dimensional matrix y (the derivative with respect
%          to BM). Typical calling sequence: y=s(t,x), returning
%          the value of s at time t in position x.
%          s can be also defined to work on multiple values of x - it is
%          much quicker then for multiple initial conditions and in that
%          case s should return a 2-D array, where the n x n-dimensional
%          matrices are joined horizontally. s is returning n x n*n_trials
%       x0: initial value where integration starts (n-dimensional
%           vector or set of initial conditions as n_dim x n_trials
%           dimensional matrix).
%       tspan: Starting time and end time for integration.
%              Integration has to run from time t=tspan(1) to time
%              t=tspan(2).
%       N: number of steps for integration. The integration
%          stepsize h=(tspan(2)-tspan(1))/N should be small.
%       corr_mx (optional): correlation matrix between Brownian
%                           motions,
%       method (optional): differentiating method:
%                           -euler - works for everything (default)
%                           -heun works for everything (slightly slower than euler)
%                           -milstein - works only for diagonal noise!
%       For the methods with diagonal noise, what's out of diagonal will be omitted
%   Outputs:
%       xend: result of integration at t=tspan(2).
%       t: vector of times at which intermediate values have been
%          computed (this should have N + 1 entries).
%       xt: intermediate values (n_dim x (N + 1)-array).
%           xt(:,k) should be the solution at t(k) or
%           size(xt) = [n_dim,n_trials,(N+1)] for multiple initials.


n_dim = size(x0,1);
n_trials = size(x0,2);

default = {'corr_mx',eye(n_dim),'method','euler'};
options = MySetOptions(default,varargin);

h=(tspan(2)-tspan(1))/N;
t = tspan(1):h:tspan(2);
xt = NaN(size(x0,1),size(x0,2),N+1);
xt(:,:,1) = x0;
W = correlated_BM(options.corr_mx,h,n_trials,N);

% m and s given by the user might work on multiple initial conditions or only on one.
% If they work only on one, there is a need to make them work on multiple conditions:
if size(m(0,x0),2)==1 && size(x0,2)>1
    M = expand_m_in_2nd_dim(m);
else
    M = m;
end
if size(s(0,x0),3)==1 && size(x0,2)>1
    S = expand_s_in_3rd_dim(s);
else
    S = s;
end


if strcmpi(options.method,'euler')
    for i=1:N
        xt(:,:,i+1) = xt(:,:,i) + h * M(t(i),xt(:,:,i)) + multiply_S_dW( S(t(i),xt(:,:,i)) , W(:,:,i+1)-W(:,:,i) );
    end
% elseif strcmpi(options.method,'heun')
%     for i=1:N
%         x_pred = xt(:,:,i) + h * M(t(i),xt(:,:,i)) + multiply_S_dW( S(t(i),xt(:,:,i)) , W(:,:,i+1)-W(:,:,i) );
%         xt(:,:,i+1) = xt(:,:,i) + 1/2 * h * (M(t(i),xt(:,:,i)) + M(t(i)+h,x_pred)) + ...
%             + 1/2 * multiply_S_dW( S(t(i),xt(:,:,i)) + S(t(i)+h,x_pred) , W(:,:,i+1)-W(:,:,i) );
%     end
% elseif strcmpi(option.method,'runge') || strcmpi(option.method,'runge-kutta')
%     for i=1:N
%         x_pred = xt(:,:,i) + h * M(t(i),xt(:,:,i)) + h^(1/2) * diag3D();
%         xt(:,:,i+1) = xt(:,:,i) + 1/2 * h * (M(t(i),xt(:,:,i)) + M(t(i)+h,x_pred)) + ...
%             + 1/2 * multiply_S_dW( S(t(i),xt(:,:,i)) + S(t(i)+h,x_pred) , W(:,:,i+1)-W(:,:,i) );
%     end
elseif strcmpi(options.method,'milstein')
    for i=1:N
        for j=1:n_trials
            xt(:,j,i+1) = xt(:,j,i) + h * M(t(i),xt(:,j,i)) +  S(t(i),xt(:,j,i)) * (W(:,j,i+1)-W(:,j,i))  + ...
                + 1/2 * diag(S(t(i),xt(:,j,i))) .* diag(MyJacobian(@(x) diag(S(t(i),x)),xt(:,j,i))) .* ((W(:,j,i+1)-W(:,j,i)).^2-h);
        end
    end
end

xend = xt(:,:,end);

if size(x0,2) == 1
    xt = reshape(xt,size(x0,1),N+1);
    xend = xt(:,end);
end

end


function M = expand_m_in_2nd_dim(m)
    M = @(t,x) m_expand(t,x,m);
end
function y = m_expand(t,x,m)
    N = size(x);
    y = NaN(size(x,1),size(x,2));
    for i=1:N(2)
        y(:,i) = m(t,x(:,i));
    end
end

function S = expand_s_in_3rd_dim(s)
    S = @(t,x) s_expand(t,x,s);
end
function y = s_expand(t,x,s)
    N = size(x,2);
    y0 = s(t,x(:,1));
    y = cat(3,y0,NaN(size(y0,1),size(y0,2),N-1));
    for i=2:N
        y(:,:,i) = s(t,x(:,i));
    end
end
    
function A = multiply_S_dW(S,dW)
% The function multiplies 3-dimensional S with 2-dimensional dW in a way,
% that each square of s is multiplied by each column of dW.
% The following sizes are necessary:
% size(S) = [size(x0,1) size(x0,1) size(x0,2)]
% size(dW) = [size(x0,1) size(x0,2)]

    size_dW1 = size(dW,1);
    size_dW2 = size(dW,2);
    A = zeros(size_dW1,size_dW2);
    for i=1:size_dW2
        A(:,i) = S(:,:,i)*dW(:,i);    
    end
end