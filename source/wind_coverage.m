%% Description

% Recovers some city emission using different footprints with different
% wind coverages, in order to demonstrate the significant influence of the
% wind coverage onto sparse reconstruction.

%% Load gurobi as cvx solver
cvx_solver Gurobi_2

%% functions
get2 = @(x, i) reshape(x(i, :), size(x, 2), 1);
get3 = @(x, i) reshape(x(i, :, :), size(x, 2), size(x, 3));
get2of3 = @(x, i, j) reshape(x(i, j, :), size(x, 3));
normalize_columns = @(X)X/(diag(sqrt(diag(X'*X))));
get_normalization_vector = @(X)eye(size(X, 2))/diag(sqrt(diag(X'*X)));

%% Input data
output_path_folder = '../output/wind_coverage';

sensing_matrix_path = "../data/footprint/wind_variation.nc";

sensing_threshhold = 1e-9; % threshhold for the sensing matrix

measurement_station_names = ["A";"B"; "C"; "D"; "E"; "F"; "G"; "H";"I"; "J";"K";"L";"M";"N";"O";"P";"Q";"R";"S";"T";"U";
                "V";"W";"X";"Y";"Z";"AA";"AB";"AC";"AD";"AE";"AF";"AG";"AH";"AI";"AJ";"AK";"AL";"AI";"AJ";"AK";"AL"];

dwt_level = 3;
% amount of measurement stations per degree of freedom
percent_measurement_stations_to_use = 0.015;

start_by = 4;

use_pseudo_emission_map = false;
%% for real emission map
if use_pseudo_emission_map == false
map = "co2_munich.nc";
city_name = "Munich";

ncfile = "../data/emission_map/" + map;
co2 = ncread(ncfile, 'CO2Diffuse');
co2 = 1e-4*co2; % convert to ymol/m^2/s

longitudes = ncread(ncfile, 'longitude');
latitudes = ncread(ncfile, 'latitude');
end

%% for pseudo emission map
if use_pseudo_emission_map == true
city_name = "Pseudo";
co2 = generate_pseudo_emissions(32,32,0.8);
warning("A pseudo emission inventory is used!");
end
emission_map = co2;
            
%% Preprocessing
angle_coverage = [];
angle_names = strings(0);

nc_file = ncinfo(sensing_matrix_path);
vars = getfield(nc_file, 'Variables');

for k=1:numel(vars)
    name = getfield(vars(k), 'Name');
    letter = extractBefore(name, '_');
    if letter == 'A' & length(letter) == 1
        angle_na = extractAfter(name, '_');
        angle = str2double(angle_na);
        angle_names = [angle_names; angle_na];
        angle_coverage = [angle_coverage; angle];
    end
end


%% Data to collect


least_square_error = [];

size_x = size(emission_map, 1);
size_y = size(emission_map, 2);
emission_map_vec = reshape(emission_map, size_x * size_y, 1);


amount_measurement_stations = floor(percent_measurement_stations_to_use * size_x * size_y);
if amount_measurement_stations > length(measurement_station_names)
    error("Not enough stations available.")
end

approx_error = [];
approx_error_dwt = [];
approx_error_l2 = [];

%% Loop over all wind angles
for i=start_by:length(angle_coverage)
    angle = angle_coverage(i);
    angle_name = angle_names(i);
    
    tmp_ms_names = [];
    for j=1:amount_measurement_stations
        tmp_ms_names = [tmp_ms_names; append(measurement_station_names(j), '_', angle_name)];
    end
    collection_sensing_matrix = load_collection_footprint(sensing_matrix_path, tmp_ms_names, size_x, size_y);


    % create the sensing matrix
    sensing_matrix = [];
    for j=1:amount_measurement_stations
        tmp = get3(collection_sensing_matrix, j);
        sensing_matrix = [sensing_matrix; tmp];
    end
    amount_measurements = size(sensing_matrix, 1);
    sensing_matrix(sensing_matrix < 0) = 0;
    

    T = threshold_sensing_matrix(sensing_matrix,sensing_threshhold);
    sensing_matrix_tilde = sensing_matrix * transpose(T);
    

    obs_vec = double(sensing_matrix*emission_map_vec);
    %obs_vec_2 = double(sensing_matrix_tilde*T*emission_map_vec);
    
    max_disturbtion = sensing_threshhold*max(emission_map_vec);
    %% reconstruction L1
    
    emission_L1_vec = optimizeL1(sensing_matrix_tilde, obs_vec);
    if any(isnan(emission_L1_vec))
        error("nan");
        %emission_L1_vec=optimizeL1_noise(sensing_matrix_tilde, obs_vec, 1e-7);
    end
    
    emission_L1_vec = transpose(T) * emission_L1_vec;
    
    emission_L1 = reshape(emission_L1_vec, size_x, size_y);
    
    
    err = norm(emission_L1_vec - emission_map_vec, 2) / norm(emission_map_vec, 2)
    approx_error = [approx_error; err];
    
     %% Wavelet sensing matrix
    [dwt_sensing_mat, S, wave_size] = DWT_matrix_rows(sensing_matrix, size_x, size_y, dwt_level);

    T_dwt = threshold_sensing_matrix(dwt_sensing_mat, sensing_threshhold);

    sensing_matrix_tilde_dwt = dwt_sensing_mat * transpose(T_dwt);
    
    %% Wavelet CS
    emission_L1_dwt_vec = optimizeL1(sensing_matrix_tilde_dwt, obs_vec);
    if any(isnan(emission_L1_dwt_vec))
        error("nan");
        emission_L1_dwt_vec = optimizeL1_noise(sensing_matrix_tilde_dwt, obs_vec, 1e-7);
   end
    emission_L1_dwt_vec = transpose(T_dwt) * emission_L1_dwt_vec;
    
    emission_L1_by_DWT = waverec2(emission_L1_dwt_vec, S, 'haar');

    emission_L1_vec_by_DWT = reshape(emission_L1_by_DWT, size_x * size_y, 1);
    
    err = norm(emission_L1_vec_by_DWT - emission_map_vec, 2) / norm(emission_map_vec, 2);
    approx_error_dwt = [approx_error_dwt; err];
    
    %% reconstruction L2
    emission_l2_vec = pinv(sensing_matrix_tilde, 1e-8) * obs_vec;
    emission_l2_vec = transpose(T) * emission_l2_vec;
    
    emission_l2 = reshape(emission_l2_vec, size_x, size_y);
    
    % imagesc(emission_l2);
    
    err = norm(emission_l2_vec - emission_map_vec, 2) / norm(emission_map_vec, 2)
    approx_error_l2 = [approx_error_l2; err];
    
end


angle_coverage_degree = 360 * angle_coverage(start_by:end)/(2 * pi);
%% Generate plotting path
output_path = strcat(output_path_folder, "/", city_name, "/", "s_", string(size_x * size_y), "ms_", string(percent_measurement_stations_to_use), "/");

mkdir(output_path);


%% Plotting
h = figure();
hold on
plot(angle_coverage_degree, approx_error, 'x--', 'linewidth', 1.5)
plot(angle_coverage_degree, approx_error_dwt, 'x--', 'linewidth', 1.5)
plot(angle_coverage_degree, approx_error_l2, 'x--', 'linewidth', 1.5)
lgd = legend("SR", "SR DWT", "LS");
lgd.Location = 'northwest';
xlabel("wind coverage [°]",'FontWeight','bold')
ylabel("Rel. l_2 error",'FontWeight','bold')
title(city_name);
xlim([80 360])
%set(gca,'YScale', 'log');
ylim([0 0.8])
grid minor
set(findall(gcf,'-property','FontSize'),'FontSize',12)
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,strcat(output_path, 'wind_variation'),'-dpdf','-r0')
hold off