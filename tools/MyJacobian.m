function df = MyJacobian(f,x,h)
%% function df = MyJacobian(f,x,h)
%
%   Inputs
%       f: function to be differentiated,
%       x: point where Jacobian is taken,
%       h (optional): parameter for finite differences [1e-6].
%   Output
%       df: m by n matrix, Jacobian of f in x.
    if nargin == 2
        h = 1e-6;
    end
    
    H = h*eye(length(x));
    di = @(F,x,H,i) (F(x+H(:,i))-F(x-H(:,i))) / (2*H(1,1));
    df = zeros(length(f(x)),length(x));
    for k = 1:length(x)
        df(:,k) = di(f,x,H,k);
    end
end