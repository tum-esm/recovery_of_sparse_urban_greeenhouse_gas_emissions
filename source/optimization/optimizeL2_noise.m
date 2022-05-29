function x = optimizeL2_noise(A,b,noise_level)
%OPTIMIZEL0 Summary of this function goes here
%   Detailed explanation goes here
n = size(A, 2);
cvx_begin;
    cvx_precision best
    variable x_l2(n);
    minimize( norm(x_l2, 2) );
    subject to
    0.999 * sum_square(A * x_l2 - b) <= noise_level;
cvx_end;
x = x_l2;
end

