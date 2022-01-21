function [collection_sensing_matrix] = load_collection_footprint(sensing_matrix_path,measurement_station_names, size_x, size_y)
%LOAD_COLLECTION_FOOTPRINT Summary of this function goes here
%   Detailed explanation goes here
amount_measurement_stations = size(measurement_station_names, 1);
tmp = ncread(sensing_matrix_path, measurement_station_names(1));
foot_size_x = size(tmp, 2);
foot_size_y = size(tmp, 3);

collection_sensing_matrix = [];
% reshape footprints to the size of the CO2 map
for i=1:amount_measurement_stations
    tmp = ncread(sensing_matrix_path, measurement_station_names(i));
    tmp = reshape(tmp, [size(tmp, 1), foot_size_x * foot_size_y]);
    tmp = reshape_footprint(tmp, foot_size_x, foot_size_y, size_x, size_y);
    tmp = reshape(tmp, [1, size(tmp, 1), size_x * size_y]);
    tmp(tmp < 0) = 0;
    collection_sensing_matrix = [collection_sensing_matrix; tmp];
end
end

