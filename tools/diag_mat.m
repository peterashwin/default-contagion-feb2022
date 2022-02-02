function y = diag_mat(A)
%UNTITLED2 make diagonal matrices from a matrix column by column
%   A - matrix of any size
%   y - stacks matrices which are diagonals of each column of A horizontally
d1 = size(A,1);
d2 = size(A,2);
y = zeros(d1,d1*d2);
for d=1:d2
    y(:,d1*(d-1)+1:d1*d) = diag(A(:,d));
end
end

