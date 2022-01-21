function [new_footprint] = reshape_footprint(footprint, old_size_x, old_size_y, new_size_x, new_size_y)
%RESHAPE_FOOTPRINT Summary of this function goes here
%   Detailed explanation goes here
scale_x = old_size_x / new_size_x;
scale_y = old_size_y / new_size_y;
correction_factor = scale_x * scale_y;
fprintf("correction factor for footprint is %f\n", correction_factor);
amount_measurements = size(footprint, 1);
new_footprint = zeros(amount_measurements, new_size_x * new_size_y);
for i=1:amount_measurements
   foot = footprint(i,:);
   foot = reshape(foot, old_size_x, old_size_y);
   foot = correction_factor * imresize(foot, [new_size_x new_size_y]);
   new_footprint(i,:) = reshape(foot, 1, new_size_x * new_size_y);
end
end

