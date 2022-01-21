function x = optimizeL1_noise(A,b, noise_level)
%OPTIMIZEL0 Summary of this function goes here
%   Detailed explanation goes here
n = size(A, 2);
cvx_begin;
    cvx_precision best
    variable x_l1(n);
    minimize( norm(x_l1, 1) );
    subject to
    0.999 * sum_square(A * x_l1 - b) <= noise_level; %% for numerical reasons 0.999
cvx_end;
x = x_l1;
end

