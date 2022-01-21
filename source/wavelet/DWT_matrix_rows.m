function [new_mat,S,wave_size] = DWT_matrix_rows(mat, size_x, size_y, level)
%DWT_MATRIX_ROWS Summary of this function goes here
%   Detailed explanation goes here
m = size(mat, 1);
n = size(mat, 2);
tmp = zeros(size_x, size_y);
tmp(1, 1) = 1;
[tmp, S] = wavedec2(tmp, level, 'haar');
wave_size = size(tmp, 2);
disp(wave_size);
new_mat = zeros(m, wave_size);
for i=1:m
    [new_mat(i, :), S] = wavedec2(reshape(mat(i,:), size_x, size_y), level, 'haar');
end
end

