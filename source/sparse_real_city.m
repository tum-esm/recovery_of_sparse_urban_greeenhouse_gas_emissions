%% Description

% Comparison of spatial vs wavelet domain for reconstructing emissions
% depending on their sparseness. To get emissions close to real world
% emissions, a city inventory is used which is then sparsified.

%% Load gurobi as cvx solver
cvx_solver Gurobi_2

%% functions
get2 = @(x, i) reshape(x(i, :), size(x, 2), 1);
get3 = @(x, i) reshape(x(i, :, :), size(x, 2), size(x, 3));
get2of3 = @(x, i, j) reshape(x(i, j, :), size(x, 3));
normalize_columns = @(X)X/(diag(sqrt(diag(X'*X))));
get_normalization_vector = @(X)eye(size(X, 2))/diag(sqrt(diag(X'*X)));


%% Input data
output_path_folder = '../output/sparse_real_city/';

sensing_matrix_path = "../data/footprint/synthtetic.nc";
sensing_threshhold = 1e-9; % threshhold for the sensing matrix

measurement_station_names = ["A"; "B"; "C"; "D"; "E"; "F"; "G"; "H"; 
                "I"; "J";"K";"L";"M";"N";"O";"P";"Q";"R";"S";"T";"U";
                "V";"X";"Y";"Z";"AA";"AB";"AC";"AD";"AE";"AF";"AG";"AH";"AI";
                "AJ";"AK";"AL";"AM";"AN"];
amount_measurement_stations = size(measurement_station_names, 1);

dwt_level = 3;

% amount of measurement stations per degree of freedom
percent_measurement_stations_to_use = 0.015;

use_pseudo_emission_map = true;
%% for real emission map
if use_pseudo_emission_map == false
map = "co2_munich.nc";
city_name = "Munich";

ncfile = "../data/emission_map/" + map;
co2 = ncread(ncfile, 'CO2Diffuse');
co2 = 1e-4 * co2; % convert to ymol/m^2/s

longitudes = ncread(ncfile, 'longitude');
latitudes = ncread(ncfile, 'latitude');
end

%% for pseudo emission map
if use_pseudo_emission_map == true
city_name = "Pseudo";
co2 = generate_pseudo_emissions(32,32,0.8);
warning("A pseudo emission inventory is used!");
end


%% Data to collect
s_vec = [];
s_of_DoF_vec = [];
error = [];
error_spatial = [];
error_L2 = [];
noise_L2 = [];

%% Preprocessing

size_x = size(co2, 1);
size_y = size(co2, 2);
co2_vec = reshape(co2, size_x * size_y, 1);



collection_sensing_matrix = load_collection_footprint(sensing_matrix_path, measurement_station_names, size_x, size_y);

amount_measurement_stations_to_use = floor(size_x * size_y * percent_measurement_stations_to_use)
if amount_measurement_stations_to_use > amount_measurement_stations
    error('not enough measurement stations in the footprint!');
end

fprintf("Amount measurement stations used is %d\n", amount_measurement_stations_to_use);


% create the sensing matrix
sensing_matrix = [];
for j=1:amount_measurement_stations_to_use
    tmp = get3(collection_sensing_matrix, j);
    sensing_matrix = [sensing_matrix; tmp];
end
amount_measurements = size(sensing_matrix, 1);
sensing_matrix(sensing_matrix < 0) = 0;

%% Create approx sensing matrix
T = threshold_sensing_matrix(sensing_matrix,sensing_threshhold);

sensing_matrix_tilde = sensing_matrix * transpose(T);

%% Wavelet sensing matrix
[dwt_sensing_mat, S, wave_size] = DWT_matrix_rows(sensing_matrix, size_x, size_y, dwt_level);

T_dwt = threshold_sensing_matrix(dwt_sensing_mat, sensing_threshhold);

sensing_matrix_tilde_dwt = dwt_sensing_mat * transpose(T_dwt);


%% Loop over s
for s_perc=0.01:0.045:1
    s = floor(s_perc * size_x * size_y);
    
    emission_map = reshape(get_LargestValueAbsolute(co2_vec, s), size_x, size_y);
    emission_map_vec = reshape(emission_map, size_x * size_y, 1);
    
    real_s = sum(emission_map_vec ~= 0);

    s_vec = [s_vec; real_s];
    s_of_DoF_vec = [s_of_DoF_vec; real_s/(size_x * size_y)];
    
%% Create observation
    obs_vec = double(sensing_matrix * emission_map_vec);
    
%% CS reconstruction
    emission_L1_vec = optimizeL1(sensing_matrix_tilde, obs_vec);

    emission_L1_vec = transpose(T) * emission_L1_vec;


    emission_L1 = reshape(emission_L1_vec, size_x, size_y);

    % imagesc(emission_L1);
    err_spatial = norm(emission_L1_vec - emission_map_vec, 2) / norm(emission_map_vec, 2)

    error_spatial = [error_spatial; err_spatial];

%% Wavelet CS
    emission_L1_dwt_vec = optimizeL1(sensing_matrix_tilde_dwt, obs_vec);

    emission_L1_dwt_vec = transpose(T_dwt) * emission_L1_dwt_vec;

    % We recover emissions in the wavelet domain directly
    emission_L1 = waverec2(emission_L1_dwt_vec, S, 'haar');
    % emission_L1_vec = emission_L1_dwt_vec;

    emission_L1_vec = reshape(emission_L1, size_x * size_y, 1);


    err = norm(emission_L1_vec - emission_map_vec, 2) / norm(emission_map_vec, 2)

    error = [error; err];

%% Least Square
    emission_l2_vec = pinv(sensing_matrix_tilde) * obs_vec;
    emission_l2_vec = transpose(T) * emission_l2_vec;

    emission_l2 = reshape(emission_l2_vec, size_x, size_y);

    % imagesc(emission_l2);

    err = norm(emission_l2_vec - emission_map_vec, 2) / norm(emission_map_vec, 2)
    error_L2 = [error_L2; err];


end


[s_of_DoF_vec, idx] = sort(s_of_DoF_vec);
error = error(idx);
error_L2 = error_L2(idx);


%% Generate plotting path
output_path = strcat(output_path_folder, "/", city_name, "/", "s_", string(real_s), "ms_", string(percent_measurement_stations_to_use), "/");

mkdir(output_path);

%% Plotting
h=figure()
plot(s_of_DoF_vec, error_spatial, 'x--', 'linewidth',1.5);
hold on
plot(s_of_DoF_vec, error, 'x--', 'linewidth',1.5);
plot(s_of_DoF_vec, error_L2, 'x--', 'linewidth',1.5);
xlabel("Rel. sparsity",'FontWeight','bold');
ylabel("Rel. l_2 error",'FontWeight','bold');
lgd = legend('SR', 'SR DWT', 'LS');
lgd.Location = 'southeast';
title(city_name)
ylim([0.0 0.6]);
grid minor
set(findall(gcf,'-property','FontSize'),'FontSize',12)
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,strcat(output_path, 'sparse_real_city'),'-dpdf','-r0')

    