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
map = "co2_vienna.nc";
city_name = "Vienna";

output_path_folder = '../output/wind_coverage';

sensing_matrix_path = "../data/footprint/wind_variation.nc";

sensing_threshhold = 1e-9; % threshhold for the sensing matrix

measurement_station_names = ["A"; "B"; "C"; "D"; "E"; "F"; "G"; "H"; 
                "I"; "J";"K";"L";"M";"N";"O";"P";"Q";"R";"S";"T";"U";
                "V";"W";"X";"Y";"Z";"AA";"AB";"AC";"AD";"AE";"AF";"AG";"AH";"AI";"AJ";"AK";"AL";"AI";"AJ";"AK";"AL"];

% amount of measurement stations per degree of freedom
percent_measurement_stations_to_use = 0.015;
            
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

ncfile = "../data/emission_map/" + map;
emission_map = ncread(ncfile, 'CO2Diffuse');
emission_map = 1e-4 * emission_map;


size_x = size(emission_map, 1);
size_y = size(emission_map, 2);
emission_map_vec = reshape(emission_map, size_x * size_y, 1);


amount_measurement_stations = floor(percent_measurement_stations_to_use * size_x * size_y);
if amount_measurement_stations > length(measurement_station_names)
    error("Not enough stations available.")
end

approx_error = [];
approx_error_l2 = [];

%% Loop over all wind angles
for i=1:length(angle_coverage)
    angle = angle_coverage(i);
    angle_name = angle_names(i);
    
    sensing_matrix = [];
    for i=1:amount_measurement_stations
        tmp = ncread(sensing_matrix_path, append(measurement_station_names(i), '_', angle_name));
        tmp2 = [];
        for j=1:size(tmp, 1)
            tmp3 = reshape(tmp(j, :, :), size(tmp, 2), size(tmp, 3));
            tmp2 = [tmp2; reshape(imresize(tmp3, [size_x, size_y]), [1, size_x, size_y])];
        end
        tmp2 = reshape(tmp2, [size(tmp2, 1), size(tmp2, 2)*size(tmp2, 3)]);
        sensing_matrix = [sensing_matrix; tmp2];
    end
    sensing_matrix(sensing_matrix < 0) = 0;
    

    temp_sensing_mat = sensing_matrix;
    temp_sensing_mat(temp_sensing_mat < sensing_threshhold) = 0;
    reduction_mask = ~all(~temp_sensing_mat,1);
    T = eye(size_x * size_y);
    T = T(reduction_mask, :);
    
    sensing_matrix_tilde = sensing_matrix * transpose(T);
    

    obs_vec = double(sensing_matrix_tilde * T * emission_map_vec);
    
 
    emission_L1_vec = optimizeL1(sensing_matrix_tilde, obs_vec);
    
    emission_L1_vec = transpose(T) * emission_L1_vec;
    while any(isnan(emission_L1_vec))
       error("NAN");
   end
    
    emission_L1 = reshape(emission_L1_vec, size_x, size_y);
    
    
    err = norm(emission_L1_vec - emission_map_vec, 2) / norm(emission_map_vec, 2)
    approx_error = [approx_error; err];
    
    emission_l2_vec = pinv(sensing_matrix_tilde, 1e-8) * obs_vec;
    emission_l2_vec = transpose(T) * emission_l2_vec;
    
    emission_l2 = reshape(emission_l2_vec, size_x, size_y);
    
    % imagesc(emission_l2);
    
    err = norm(emission_l2_vec - emission_map_vec, 2) / norm(emission_map_vec, 2)
    approx_error_l2 = [approx_error_l2; err];
    
end


angle_coverage = 360 * angle_coverage/(2 * pi);
%% Generate plotting path
output_path = strcat(output_path_folder, "/", city_name, "/", "s_", string(size_x * size_y), "ms_", string(percent_measurement_stations_to_use), "/");

mkdir(output_path);


%% Plotting
h = figure();
hold on
plot(angle_coverage, approx_error, 'x--', 'linewidth', 1.5)
plot(angle_coverage, approx_error_l2, 'x--', 'linewidth', 1.5)
legend("SR", "LS")
xlabel("wind coverage [°]",'FontWeight','bold')
ylabel("Rel. l_2 error",'FontWeight','bold')
title(city_name);
xlim([0 360])
grid minor
set(findall(gcf,'-property','FontSize'),'FontSize',12)
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,strcat(output_path, 'wind_variation'),'-dpdf','-r0')
hold off