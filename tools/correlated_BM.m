function corr_W = correlated_BM(corr_mx,h,k,N)
%% corr_W = correlated_BM(corr_mx,h,n2,N)
%   Inputs:
%       corr_mx: correlation matrix (positive definite) of size n1 x n1.
%       h: size of timestep.
%       k: number of realizations.
%       N: number of steps for simulation.
%   Outputs:
%       corr_W: 3-D matrix (dim,k,N) of sets of correlated BM realizations.

d = size(corr_mx,1); %how many Brownians to correlate
C = chol(corr_mx);
corr_W = zeros(d,k,N+1);

for i=1:k
    rnd = randn(N,d);
    corr_W(:,i,:) = [zeros(1,d);h.^(1/2).*cumsum(rnd*C)]';
end

% for i=1:min([n2 5])
%     figure
%     x = reshape(corr_W(:,i,:),n1,N+1);
%     plot(0:N,x);
% end
end