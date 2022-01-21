function x = optimizeL1(A,b)
%OPTIMIZEL1 Summary of this function goes here
%   Detailed explanation goes here
n = size(A, 2);
cvx_begin;
    cvx_precision best
    variable x_l1(n);
    minimize( norm(x_l1, 1) );
    subject to
    A * x_l1 == b;
cvx_end;
x = x_l1;
end

