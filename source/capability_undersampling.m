%% Description

% Demonstrates the undersampling capability of sparse reconstruction in the
% noiseless case.

%% Load gurobi as cvx solver
cvx_solver Gurobi

%% functions
get2 = @(x, i) reshape(x(i, :), size(x, 2), 1);
get3 = @(x, i) reshape(x(i, :, :), size(x, 2), size(x, 3));
get2of3 = @(x, i, j) reshape(x(i, j, :), size(x, 3));
normalize_columns = @(X)X/(diag(sqrt(diag(X'*X))));
get_normalization_vector = @(X)eye(size(X, 2))/diag(sqrt(diag(X'*X)));


%% Parameters
maps = ["co2_berlin.nc"; 
        "co2_hamburg.nc"; 
        "co2_london.nc";
        "co2_munich_all_sources.nc";
        "co2_paris.nc";
        "co2_vienna.nc"];

map = "co2_paris.nc";
city_name = "Paris";

sensing_matrix_path = "../data/footprint/synthtetic.nc";
sensing_threshhold = 1e-9; % threshhold for the sensing matrix

measurement_station_names = ["A"; "B"; "C"; "D"; "E"; "F"; "G"; "H"; 
                "I"; "J";"K";"L";"M";"N";"O";"P";"Q";"R";"S";"T";"U";
                "V"; "W"; "X";"Y";"Z";"AA";"AB";"AC";"AD";"AE";"AF";"AG";"AH";"AI";
                "AJ";"AK";"AL";"AM";"AN";"AO";"AP";"AQ";"AR";"AS";"AT";"AU";
                "AV";"AW";"AX";"AY";"AZ";"BA";"BB";"BC";"BD";"BE";"BF";"BG";
                "BH";"BI";"BJ";"BK";"BL";"BM";"BN";"BO";"BP";"BQ";"BR";"BS";
                "BT";"BU";"BV";"BW";"BX";"BY";"BZ";"CA"];
amount_measurement_stations = size(measurement_station_names, 1);

dwt_level = 3;

% amount of measurement stations per degree of freedom
percent_measurement_stations_to_use = 0.015;

%% Data to collect
relative_errors = [];
relative_errors_normalized = [];
relative_errors_wavelet = [];

L1_relative_errors = [];
L1_relative_errors_wavelet = [];

percent_measurements = [];

L1_minimum_error = [];
L1_minimum_error_wavelet = [];

least_square_error = [];

ncfile = "../data/emission_map/" + map;
co2 = ncread(ncfile, 'CO2Diffuse');

co2 = 1e-4 * co2; % convert co2 scale

size_x = size(co2, 1);
size_y = size(co2, 2);

co2_vec = reshape(co2, size_x * size_y, 1);

collection_sensing_matrix = load_collection_footprint(sensing_matrix_path, measurement_station_names, size_x, size_y);

use_every_nth_measurement = 1;
amount_measurements_per_station = size(collection_sensing_matrix, 2);
amount_measurements_per_step = 2;
amount_measurement_stations_to_use = floor(size_x * size_y * percent_measurement_stations_to_use);
if amount_measurement_stations_to_use > amount_measurement_stations
    error('not enough measurement stations in the footprint!');
end


for i=1:amount_measurements_per_step:amount_measurements_per_station
    sensing_matrix = [];
    for j=1:amount_measurement_stations_to_use
        tmp = get3(collection_sensing_matrix, j);
        tmp = tmp(1:use_every_nth_measurement:i, :);
        sensing_matrix = [sensing_matrix; tmp];
    end
    sensing_matrix(sensing_matrix < 0) = 0;
    amount_measurements = size(sensing_matrix, 1);

    temp_sensing_mat = sensing_matrix;
    temp_sensing_mat(temp_sensing_mat < sensing_threshhold) = 0;
    reduction_mask = ~all(~temp_sensing_mat,1);
    T = eye(size_x * size_y);
    T = T(reduction_mask, :);
    
    P = transpose(T) * T; % Projector on emitters we want to dertermine
    sensing_matrix_tilde = sensing_matrix * transpose(T);

    percent_measurements = [percent_measurements; amount_measurements/(size(sensing_matrix, 2))];

    %% Delete what is not measurable on the CO2 vector
    co2_vec_tilde = co2_vec;
    
    %% Wavelet sensing matrix
    [dwt_sensing_mat, S, wave_size] = DWT_matrix_rows(sensing_matrix, size_x, size_y, dwt_level);
    
    temp_sensing_mat = dwt_sensing_mat;
    temp_sensing_mat(temp_sensing_mat < sensing_threshhold) = 0;
    dwt_reduction_mask = ~all(~temp_sensing_mat,1);
    
    T_dwt = eye(wave_size);
    T_dwt = T_dwt(dwt_reduction_mask, :);
    
    P_dwt = transpose(T_dwt) * T_dwt; % Projector on emiiters we want to dertermine
    sensing_matrix_tilde_dwt = dwt_sensing_mat * transpose(T_dwt);
    
    %% Create observation
    obs_vec = double(sensing_matrix * co2_vec);
    %% reconstruction without DWT
    
    co2_l1 = optimizeL1(sensing_matrix_tilde, obs_vec);
    
    co2_l1 = transpose(T) * co2_l1;
    
    rel_error = norm(co2_vec_tilde - co2_l1, 2)/norm(co2_vec_tilde, 2);
    relative_errors = [relative_errors; rel_error];
    
    %% Wavelet reconstruction
    co2_l1_dwt = optimizeL1(sensing_matrix_tilde_dwt, obs_vec);
    co2_l1_dwt = transpose(T_dwt) * co2_l1_dwt;
    
    rec_co2_l1_dwt = waverec2(co2_l1_dwt, S, 'haar');
    
    rec_co2_l1_dwt = reshape(rec_co2_l1_dwt, size_x * size_y, 1);
    
    rel_error_dwt = norm(co2_vec_tilde - rec_co2_l1_dwt, 2)/norm(co2_vec_tilde, 2);
    relative_errors_wavelet = [relative_errors_wavelet; rel_error_dwt];
    
    
    %% L2 Reconstruction
    co2_l2 = pinv(sensing_matrix_tilde) * obs_vec;
    co2_l2 = transpose(T) * co2_l2;
    
    ls_error = norm(co2_vec_tilde - co2_l2, 2)/norm(co2_vec_tilde, 2);
    least_square_error = [least_square_error; ls_error];
end

%% Plotting
h = figure();
hold on
title(city_name);
xlabel("m/n",'FontWeight','bold');
ylabel("Rel. error",'FontWeight','bold');
plot(percent_measurements, relative_errors, 'x--', 'linewidth',1.5);
% plot(percent_measurements, relative_errors_normalized, 'x--');
plot(percent_measurements, relative_errors_wavelet, 'x--', 'linewidth',1.5);
plot(percent_measurements, least_square_error, 'x--', 'linewidth',1.5)
legend('SR', ...
    'SR DWT',...
    "LS");
grid minor
set(findall(gcf,'-property','FontSize'),'FontSize',12)
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h, sprintf('../output/undersampling_%s', city_name),'-dpdf','-r0')
hold off

