function T = threshold_sensing_matrix(sensing_matrix,threshold)
%THRESHOLD_SENSING_MATRIX Summary of this function goes here
%   Detailed explanation goes here
T = eye(size(sensing_matrix,2));
T = T(vecnorm(sensing_matrix, 1) > threshold,:);
end

