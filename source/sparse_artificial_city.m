%% Description

% This demonstrates the exact recovery of sparse reconstruction for exactly
% sparse emissions using some footprints, even though the footprints do not
% fulfill all requirements for comrpessed sensing (run
% "analyze_sensing_matrix.m" for that).

%% functions
get2 = @(x, i) reshape(x(i, :), size(x, 2), 1);
get3 = @(x, i) reshape(x(i, :, :), size(x, 2), size(x, 3));
get2of3 = @(x, i, j) reshape(x(i, j, :), size(x, 3));
normalize_columns = @(X)X/(diag(sqrt(diag(X'*X))));
get_normalization_vector = @(X)eye(size(X, 2))/diag(sqrt(diag(X'*X)));


%% Input data
output_path_folder = '../output/sparse/';
sensing_matrix_path = "../data/footprint/synthtetic.nc";

sensing_threshhold = 1e-9; % threshhold for the sensing matrix

measurement_station_names = ["A"; "B"; "C"; "D"; "E"; "F"; "G"; "H"; 
                "I"; "J";"K";"L";"M";"N";"O";"P";"Q";"R";"S";"T";"U";
                "V";"X";"Y";"Z";"AA";"AB";"AC";"AD";"AE";"AF";"AG"];
amount_measurement_stations = size(measurement_station_names, 1);

% sparsity
s = 475;

% amount of measurement stations per degree of freedom
percent_measurement_stations_to_use = 0.015;

% size
size_x = 32;
size_y = 32;


%% Create emission map

emission_map = create_s_sparse_emission_map(s, size_x, size_y, true);

emission_map_vec = reshape(emission_map, size_x * size_y, 1);

real_s = sum(emission_map_vec ~= 0);

% Determine the amount of measurement matrices to use
amount_measurement_stations_to_use = floor(size_x * size_y * percent_measurement_stations_to_use);

collection_sensing_matrix = [];
% reshape footprints to the size of the CO2 map
for i=1:amount_measurement_stations
    tmp = ncread(sensing_matrix_path, measurement_station_names(i));
    tmp2 = [];
    for j=1:size(tmp, 1)
        tmp3 = reshape(tmp(j, :, :), size(tmp, 2), size(tmp, 3));
        tmp2 = [tmp2; reshape(imresize(tmp3, [size_x, size_y]), [1, size_x, size_y])];
    end
    tmp2 = reshape(tmp2, [1, size(tmp2, 1), size(tmp2, 2)*size(tmp2, 3)]);
    collection_sensing_matrix = [collection_sensing_matrix; tmp2];
end

amount_stations_available = size(collection_sensing_matrix, 1);

if amount_stations_available < amount_measurement_stations_to_use
    amount_measurement_stations_to_use = amount_stations_available;
end

fprintf("Amount of measurement stations used is  %d\n", amount_measurement_stations_to_use);


% create the sensing matrix
sensing_matrix = [];
for j=1:amount_measurement_stations_to_use
    tmp = get3(collection_sensing_matrix, j);
    sensing_matrix = [sensing_matrix; tmp];
end
% sensing_matrix = normalize_columns(sensing_matrix);
sensing_matrix(sensing_matrix < 0) = 0;
amount_measurements = size(sensing_matrix, 1);


% Apply threshholding to the matrix
temp_sensing_mat = sensing_matrix;
temp_sensing_mat(temp_sensing_mat < sensing_threshhold) = 0;
reduction_mask = ~all(~temp_sensing_mat,1);
T = eye(size_x * size_y);
T = T(reduction_mask, :);

sensing_matrix_tilde = sensing_matrix * transpose(T);

%% Create observations
obs_vec = double(sensing_matrix_tilde * T * emission_map_vec);


%% recosntruct L1
emission_L1_vec = optimizeL1(sensing_matrix_tilde, obs_vec);

emission_L1_vec = transpose(T) * emission_L1_vec;

emission_L1 = reshape(emission_L1_vec, size_x, size_y);


err = norm(emission_L1_vec - emission_map_vec, 2) / norm(emission_map_vec, 2);



%% reconstruction L2
emission_l2_vec = pinv(sensing_matrix_tilde, 1e-8) * obs_vec;
emission_l2_vec = transpose(T) * emission_l2_vec;

emission_l2 = reshape(emission_l2_vec, size_x, size_y);

err_l2 = norm(emission_l2_vec - emission_map_vec, 2) / norm(emission_map_vec, 2);


%% Generate output path
output_path = strcat(output_path_folder, "/s_", string(real_s), "ms_", string(percent_measurement_stations_to_use), "/");

mkdir(output_path);

%% Open file to write information to
fileID = fopen(strcat(output_path, 'info.txt'),'w');
%% Printing information
fprintf(fileID, "Info file for sparse reconstruction\n");
fprintf(fileID, "\n");
fprintf(fileID, "n is %d.\n", size_x * size_y);
fprintf(fileID, "m is %d.\n", size(sensing_matrix_tilde, 1));
fprintf(fileID, "s is %d\n", real_s);
fprintf(fileID, "Amount of measurement stations used is  %d\n", amount_measurement_stations_to_use);
fprintf(fileID, "\n");
fprintf(fileID, "The l2 reconstruction error using sparse reconstruction is %f\n", norm(emission_L1_vec - emission_map_vec, 2)/norm(emission_map_vec,2));
fprintf(fileID, "The l1 reconstruction error using sparse reconstruction is %f\n", norm(emission_L1_vec - emission_map_vec, 1)/norm(emission_map_vec,1));
fprintf(fileID, "\n");
fprintf(fileID, "The l2 reconstruction error using least squares is %f\n", norm(emission_l2_vec - emission_map_vec, 2)/norm(emission_map_vec,2));
fprintf(fileID, "The l1 reconstruction error using least squares is %f\n", norm(emission_l2_vec - emission_map_vec, 1)/norm(emission_map_vec,1));
fprintf(fileID, "\n");
%% Plotting l1 reconstructed emission map vs real emission map
max_color = max(max(max(emission_map)), max(max(emission_L1)));
h=figure('Position', [100, 100, 1024, 512]);
h1 = subplot(1,2,1);
imagesc(emission_L1, [0, max_color]);
title("reconstructed emission");
hold on
h2 = subplot(1,2,2);
imagesc(emission_map, [0, max_color]);
title("real emission");
sgtitle("Comparison L1 reconstruction to Real Emissions")
%grid minor
set(findall(gcf,'-property','FontSize'),'FontSize',12)
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,strcat(output_path, 'maps_l1_real_comparision'),'-dpdf','-r0')


%% Plotting not discovered elements vs false discovered elements

h = figure();
discovered_l1 = emission_L1;
discovered_l1(abs(discovered_l1) > 1e-1) = 1;
discovered_l1(abs(discovered_l1) < 2e-1) = 0;

discovered_real = emission_map;
discovered_real(abs(discovered_real) > 1e-2) = 1;
discovered_real(abs(discovered_real) < 1e-1) = 0;

% 1 => not dicovered
% 0 => discovered
% -1 => false discovered
false_discoveries_l1 = discovered_real - discovered_l1;

imagesc(false_discoveries_l1, [-1, 1]);
title("False Discovery Map");
hold on

xlabel("west - east [km]",'FontWeight','bold');
ylabel("south - north [km]",'FontWeight','bold');
grid minor
set(findall(gcf,'-property','FontSize'),'FontSize',12)
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,strcat(output_path, 'false_discoveries_BP'),'-dpdf','-r0')

%% Plotting not discovered elements vs false discovered elements L2

h = figure();
discovered_l2 = emission_l2;
discovered_l2(abs(discovered_l2) > 1e-1) = 1;
discovered_l2(abs(discovered_l2) < 2e-1) = 0;


discovered_real = emission_map;
discovered_real(abs(discovered_real) > 1e-2) = 1;
discovered_real(abs(discovered_real) < 1e-1) = 0;

% 1 => not dicovered
% 0 => discovered
% -1 => false discovered
false_discoveries_l2 = discovered_real - discovered_l2;

imagesc(false_discoveries_l2, [-1, 1]);
title("False Discovery Map");
hold on

xlabel("west - east [km]",'FontWeight','bold');
ylabel("south - north [km]",'FontWeight','bold');
grid minor
set(findall(gcf,'-property','FontSize'),'FontSize',12)
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,strcat(output_path, 'false_discoveries_least_square'),'-dpdf','-r0')

%% Printing more information
% Emissions falsely recovered
emission_strength_false_discovered_l1 = norm(reshape(emission_L1(false_discoveries_l1 < 0), [],1), 1) / norm(emission_L1_vec, 1);
emission_strength_false_discovered_l2 = norm(reshape(emission_l2(false_discoveries_l2 < 0), [],1), 1)/ norm(emission_l2_vec, 1);

fprintf(fileID, "Fraction of false discovered emissions for 11 reconstruction divided by total emissions is %f\n", emission_strength_false_discovered_l1);
fprintf(fileID, "Fraction of false discovered emissions for 12 reconstruction divided by total emissions is %f\n", emission_strength_false_discovered_l2);
fprintf(fileID, "\n");
fprintf(fileID, "Amount false discovered emitters for l1 reconstruction %f\n", sum(sum(false_discoveries_l1 < 0)));
fprintf(fileID, "Amount false discovered emitters for l2 reconstruction %f\n", sum(sum(false_discoveries_l2 < 0)));
fprintf(fileID, "\n");
fclose(fileID);