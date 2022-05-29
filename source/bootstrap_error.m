%% Description

% Bootstrap estimation of variance and other stuff.

%% functions
get2 = @(x, i) reshape(x(i, :), size(x, 2), 1);
get3 = @(x, i) reshape(x(i, :, :), size(x, 2), size(x, 3));
get2of3 = @(x, i, j) reshape(x(i, j, :), size(x, 3));
normalize_columns = @(X)X/(diag(sqrt(diag(X'*X))));
get_normalization_vector = @(X)eye(size(X, 2))/diag(sqrt(diag(X'*X)));
logn = @(a, n) log(a) / log(n);

%% INPUT data

% Amount of iterations
N = 250;


% SNR_list = [400, 7, 1];
% SNR_list = [7];
SNR_list = [100]; % 100 is 20 dB

amount_filtering = -1; % -1: Detect automatic how often filtering is possible

output_path_folder = '../output/case_studies';

use_pseudo_emission_map = true;
%% for real emission map
if use_pseudo_emission_map == false
map = "co2_munich.nc";
city_name = "Munich";

ncfile = "../data/emission_map/" + map;
co2 = ncread(ncfile, 'CO2Diffuse');
co2 = 1e-4 * co2; % convert units
end

%% for pseudo emission map
if use_pseudo_emission_map == true
city_name = "Pseudo";
co2 = generate_pseudo_emissions(32,32,0.8);
warning("A pseudo emission inventory is used!");
end

%% Preprocessing data
sensing_matrix_path = "../data/footprint/synthtetic.nc";
sensing_threshhold = 1e-9; % threshhold for the sensing matrix

measurement_station_names = ["A"; "B"; "C"; "D"; "E"; "F"; "G"; "H"; 
                "I"; "J";"K";"L";"M";"N";"O";"P";"Q";"R";"S";"T";"U";
                "V";"X";"Y";"Z";"AA";"AB";"AC";"AD";"AE";"AF";"AG"];
amount_measurement_stations = size(measurement_station_names, 1);

amount_SNR = length(SNR_list);

dwt_level = 3;

% amount of measurement stations per degree of freedom
percent_measurement_stations_to_use = 0.03;


size_x = size(co2, 1);
size_y = size(co2, 2);
co2_vec = reshape(co2, size_x * size_y, 1);

emission_map = co2;
emission_map_vec = co2_vec;


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

use_every_nth_measurement = 1;
amount_measurements_per_station = size(collection_sensing_matrix, 2);
amount_measurements_per_step = 4;
amount_measurement_stations_to_use = floor(size_x * size_y * percent_measurement_stations_to_use)
if amount_measurement_stations_to_use > amount_measurement_stations
    % error('not enough measurement stations in the footprint!');
    amount_measurement_stations_to_use = amount_measurement_stations;
end


% create the sensing matrix
sensing_matrix = [];
for j=1:amount_measurement_stations_to_use
    tmp = get3(collection_sensing_matrix, j);
    sensing_matrix = [sensing_matrix; tmp];
end
amount_measurements = size(sensing_matrix, 1);
sensing_matrix(sensing_matrix < 0) = 0;

%% Create approx sensing matrix
temp_sensing_mat = sensing_matrix;
temp_sensing_mat(temp_sensing_mat < sensing_threshhold) = 0;
reduction_mask = ~all(~temp_sensing_mat,1);

T = eye(size_x * size_y);
T = T(reduction_mask, :);

sensing_matrix_tilde = sensing_matrix * transpose(T);

%% Wavelet sensing matrix
[dwt_sensing_mat, S, wave_size] = DWT_matrix_rows(sensing_matrix, size_x, size_y, dwt_level);

T_dwt = eye(wave_size);

sensing_matrix_tilde_dwt = dwt_sensing_mat * transpose(T_dwt);

%% Create observation
obs_vec = double(sensing_matrix * emission_map_vec);
amount_observations = size(obs_vec, 1);

%% Generate plotting path
output_path = strcat(output_path_folder, "/", city_name, "/", "s_", ...
                     string(size_x * size_y), "ms_", string(percent_measurement_stations_to_use),...
                     "_bootstrap", "/");

mkdir(output_path);

%% Filtering functions
filter_mat = @(Mat, new_x, new_y) reshape( sum(sum(reshape(Mat, 2, new_x, 2, new_y), 1), 3), new_x, new_y);
% filter_matn = @(Mat, new_x, new_y, n) reshape( sum(sum(reshape(Mat, n, new_x, n, new_y), 1), 3), new_x, new_y);
filter_matn = @(Mat, new_x, new_y, n) smoothing(Mat, n);

x_extension = ceil(log2(size_x));
y_extension = ceil(log2(size_y));

if amount_filtering == -1
    amount_filtering = min(x_extension, y_extension);
end

new_amount_filtering = max(size_x, size_y);

%% Bootstrap integral
reconstruct_func = @(e)optimizeL1_noise(sensing_matrix_tilde, obs_vec + e, sum_square(e));
reconstruct_dwt_func = @(e)optimizeL1_noise(sensing_matrix_tilde_dwt, obs_vec + e, sum_square(e));
reconstruct_lsq_func = @(e)optimizeL2_noise(sensing_matrix_tilde, obs_vec + e, sum_square(e));

signal_power = sum_square(obs_vec) / amount_observations;

amount_captures = 10;
[map_sorted, idx_map] = sort(emission_map_vec, 'descend');
captured_points = 15.^(linspace(1, logn(size_x * size_y, 15), amount_captures));

captured_mean = zeros(amount_captures, amount_SNR);
captured_squared = zeros(amount_captures, amount_SNR);

captured_dwt_mean = zeros(amount_captures, amount_SNR);
captured_dwt_squared = zeros(amount_captures, amount_SNR);

captured_lsq_mean = zeros(amount_captures, amount_SNR);
captured_lsq_squared = zeros(amount_captures, amount_SNR);

expect_rec = zeros(size_x * size_y, amount_SNR);
expect_rec_squared = zeros(size_x * size_y, size_x * size_y, amount_SNR);

expect_rec_dwt = zeros(size_x * size_y, amount_SNR);
expect_rec_dwt_squared = zeros(size_x * size_y, size_x * size_y, amount_SNR);

new_error_filter = zeros(amount_filtering, N, amount_SNR);
new_error_filter_l1 = zeros(amount_filtering, N, amount_SNR);

new_error_filter_dwt = zeros(amount_filtering, N, amount_SNR);
new_error_filter_l1_dwt = zeros(amount_filtering, N, amount_SNR);

new_error_filter_lsq = zeros(amount_filtering, N, amount_SNR);
new_error_filter_l1_lsq = zeros(amount_filtering, N, amount_SNR);

for i=1:N
    for s=1:amount_SNR
        SNR = SNR_list(s);
        % Calculate sigma of noise
        sigma_noise = (signal_power/SNR)^0.5;
       e = sigma_noise * randn(amount_observations, 1);
       
       % reconstruction CS
       rec = transpose(T) * reconstruct_func(e);
       while any(isnan(rec))
           error("NAN");
       end
       
       % reconstruction CS in wavelet domain
       rec_dwt = transpose(T_dwt) * reconstruct_dwt_func(e);
       rec_dwt = reshape(waverec2(rec_dwt, S, 'haar'), size_x * size_y, 1);
       while any(isnan(rec_dwt))
           error("NAN");
       end
       
       % reconstruction least square
       rec_lsq = transpose(T) * reconstruct_lsq_func(e);

       rec_map = reshape(rec, size_x, size_y);
       rec_dwt_map = reshape(rec_dwt, size_x, size_y);
       rec_lsq_map = reshape(rec_lsq, size_x, size_y);

       %% Create filtered error
       filtered_map = zeros(size_x, size_y);
       filtered_dwt_map = zeros(size_x, size_y);
       filtered_lsq_map = zeros(size_x, size_y);
       filtered_emission_map = zeros(size_x, size_y);

       filtered_map(1:size_x, 1:size_y) = rec_map;
       filtered_dwt_map(1:size_x, 1:size_y) = rec_dwt_map;
       filtered_lsq_map(1:size_x, 1:size_y) = rec_lsq_map;
       filtered_emission_map(1:size_x, 1:size_y) = emission_map;
       
       err_l2 = norm(rec - emission_map_vec, 2)/norm(emission_map_vec, 2);
       err_l1 = norm(rec - emission_map_vec, 1)/norm(emission_map_vec, 1);
       
       err_l2_dwt = norm(rec_dwt - emission_map_vec, 2)/norm(emission_map_vec, 2);
       err_l1_dwt = norm(rec_dwt - emission_map_vec, 1)/norm(emission_map_vec, 1);
       
       err_l2_lsq = norm(rec_lsq - emission_map_vec, 2)/norm(emission_map_vec, 2);
       err_l1_lsq = norm(rec_lsq - emission_map_vec, 1)/norm(emission_map_vec, 1);

       compr_x = 2^x_extension;
       compr_y = 2^y_extension;
       
       % new filtering approach
       new_error_filter(1, i, s) = err_l2;
       new_error_filter_l1(1, i, s) = err_l1;
       
       new_error_filter_dwt(1, i, s) = err_l2_dwt;
       new_error_filter_l1_dwt(1, i, s) = err_l1_dwt;
       
       new_error_filter_lsq(1, i, s) = err_l2_lsq;
       new_error_filter_l1_lsq(1, i, s) = err_l1_lsq;
       
       for j=2:new_amount_filtering
            new_x_extension = ceil(size_x / j);
            new_y_extension = ceil(size_y / j);
            
            compr_x = j * new_x_extension;
            compr_y = j * new_y_extension;
            
            compr_x  = compr_x/j;
            compr_y = compr_y/j;
            
            filtered_map = zeros(size_x, size_y);
            filtered_dwt_map = zeros(size_x, size_y);
            filtered_lsq_map = zeros(size_x, size_y);
            filtered_emission_map = zeros(size_x, size_y);
            
            filtered_map(1:size_x, 1:size_y) = rec_map;
            filtered_dwt_map(1:size_x, 1:size_y) = rec_dwt_map;
            filtered_lsq_map(1:size_x, 1:size_y) = rec_lsq_map;
            filtered_emission_map(1:size_x, 1:size_y) = emission_map;
           
            filtered_map = filter_matn(filtered_map, compr_x, compr_y, j);
            filtered_dwt_map = filter_matn(filtered_dwt_map, compr_x, compr_y, j);
            filtered_lsq_map = filter_matn(filtered_lsq_map, compr_x, compr_y, j);
            filtered_emission_map = filter_matn(filtered_emission_map, compr_x, compr_y, j);
            
            filtered_map_vec = reshape(filtered_map, size_x * size_y, 1);
            filtered_map_dwt_vec = reshape(filtered_dwt_map, size_x * size_y, 1);
            filtered_map_lsq_vec = reshape(filtered_lsq_map, size_x * size_y, 1);
            filtered_emission_map_vec = reshape(filtered_emission_map, size_x * size_y, 1);
            
            new_error_filter(j,i, s) = norm(filtered_map_vec - filtered_emission_map_vec, 2)/norm(filtered_emission_map_vec, 2);
            new_error_filter_l1(j,i, s) = norm(filtered_map_vec - filtered_emission_map_vec, 1)/norm(filtered_emission_map_vec, 1);
           
            new_error_filter_dwt(j,i, s) = norm(filtered_map_dwt_vec - filtered_emission_map_vec, 2)/norm(filtered_emission_map_vec, 2);
            new_error_filter_l1_dwt(j,i, s) = norm(filtered_map_dwt_vec - filtered_emission_map_vec, 1)/norm(filtered_emission_map_vec, 1);
            
            new_error_filter_lsq(j,i, s) = norm(filtered_map_lsq_vec - filtered_emission_map_vec, 2)/norm(filtered_emission_map_vec, 2);
            new_error_filter_l1_lsq(j,i, s) = norm(filtered_map_lsq_vec - filtered_emission_map_vec, 1)/norm(filtered_emission_map_vec, 1);
       end
       
       %% Calculate captured Emissions
        [rec_sorted, idx_l1] = sort(rec, 'descend');
        [rec_dwt_sorted, idx_l1_dwt] = sort(rec_dwt, 'descend');
        [rec_lsq_sorted, idx_l1_lsq] = sort(rec_lsq, 'descend');

        for k=1:amount_captures
            j = floor(captured_points(k));
            equal_same = j;
            equal_l1 = length(intersect(idx_map(1:j), idx_l1(1:j))) / equal_same;
            equal_l1_dwt = length(intersect(idx_map(1:j), idx_l1_dwt(1:j))) / equal_same;
            equal_l1_lsq = length(intersect(idx_map(1:j), idx_l1_lsq(1:j))) / equal_same;

            captured_mean(k, s) = captured_mean(k, s) + (1/N) * equal_l1;
            captured_squared(k, s) = captured_squared(k, s) + equal_l1^2;
            
            captured_dwt_mean(k, s) = captured_dwt_mean(k, s) + (1/N) * equal_l1_dwt;
            captured_dwt_squared(k, s) = captured_dwt_squared(k, s) + equal_l1_dwt^2;
            
            captured_lsq_mean(k, s) = captured_lsq_mean(k, s) + (1/N) * equal_l1_lsq;
            captured_lsq_squared(k, s) = captured_lsq_squared(k, s) + equal_l1_lsq^2;
        end
       %% Expectation value and squared for covariance
       expect_rec(:, s) = expect_rec(:, s) + (1/N) * rec;
       expect_rec_squared(:,:,s) = expect_rec_squared(:,:,s) + (rec * rec');
       
       expect_rec_dwt(:, s) = expect_rec_dwt(:, s) + (1/N) * rec_dwt;
       expect_rec_dwt_squared(:,:,s) = expect_rec_dwt_squared(:,:,s) + (rec_dwt * rec_dwt');
    end
end

%% Calculate covariance
% cov(X) = 1/N-1 * \sum_i x_i^2 - N/N-1 E[X]^2
cov_rec = (1/(N-1)) * squeeze(expect_rec_squared(:,:,1)) - (N/(N-1)) * (squeeze(expect_rec(:,1)) * squeeze(expect_rec(:,1))');
cov_rec = squeeze(cov_rec);

var_rec = diag(cov_rec);

relative_var_rec = pinv(diag(emission_map_vec)) * var_rec;

var_map = reshape(var_rec, size_x, size_y);
relative_var_map = reshape(relative_var_rec, size_x, size_y);

%% Calculate covariance DWT
cov_rec_dwt = (1/(N-1)) * squeeze(expect_rec_dwt_squared(:,:,1)) - (N/(N-1)) * (squeeze(expect_rec_dwt(:,1)) * squeeze(expect_rec_dwt(:,1))');
cov_rec_dwt = squeeze(cov_rec_dwt);

var_rec_dwt = diag(cov_rec_dwt);

var_map_dwt = reshape(var_rec_dwt, size_x, size_y);

%% Generate correlation matrix

std_vec = sqrt(var_rec);
cor_mat = pinv(diag(std_vec)) * cov_rec * pinv(diag(std_vec));

% DWT
std_vec_dwt = sqrt(var_rec_dwt);
cor_mat_dwt = pinv(diag(std_vec_dwt)) * cov_rec * pinv(diag(std_vec_dwt));


%%  Generate distance matrix
x = 1:size_x;
y = 1:size_y;

[X,Y] = meshgrid(x,y);

dist_func = @(pos_x, pos_y) sqrt((X - pos_x).^2 + (Y - pos_y).^2);
distance_matrix = zeros(size_x * size_y, size_x * size_y);
count = 1;
for i_x = 1:size_x
    for i_y = 1:size_y
        distance_matrix(count, :) = reshape(dist_func(i_x, i_y), size_x * size_y, 1);
        count = count + 1;
    end
end


%% Store All data
save(strcat(output_path, 'project_results.mat'));
% load(strcat('../figures/case_studies/Munich/s_528ms_0.015/bootstrap/', 'project_results.mat'));

%% Start with Plotting things
%% Plot covariance
h = figure();
imagesc(sqrt(var_map));
%surf(sqrt(var_map));
set(gca,'ColorScale','log');
xlabel("west - east [km]",'FontWeight','bold');
ylabel("south - north [km]",'FontWeight','bold');
colorbar;
grid minor
set(findall(gcf,'-property','FontSize'),'FontSize',12)
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,strcat(output_path, 'variance_map'),'-dpdf','-r0')

%% Plot relative variance
h = figure();
imagesc(sqrt(relative_var_map));
set(gca,'ColorScale','log');
xlabel("west - east [km]",'FontWeight','bold');
ylabel("south - north [km]",'FontWeight','bold');
colorbar;
grid minor
set(findall(gcf,'-property','FontSize'),'FontSize',12)
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,strcat(output_path, 'rel_variance_map'),'-dpdf','-r0')

%% Plot Covariance DWT
h = figure();
imagesc(sqrt(var_map_dwt));
set(gca,'ColorScale','log');
xlabel("west - east [km]",'FontWeight','bold');
ylabel("south - north [km]",'FontWeight','bold');
colorbar;
grid minor
set(findall(gcf,'-property','FontSize'),'FontSize',12)
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,strcat(output_path, 'variance_map_dwt'),'-dpdf','-r0')

%% Plot distance vs Correlation

vectorized_distance = reshape(distance_matrix, size_x * size_y * size_x * size_y, 1);
vectorized_correlation = reshape(abs(cor_mat), size_x * size_y * size_x * size_y, 1);

h = figure();
hold on
xlabel("Distance",'FontWeight','bold');
ylabel("Correlation",'FontWeight','bold');
scatter(vectorized_distance, vectorized_correlation);
p = polyfit(vectorized_distance, vectorized_correlation, 5);
x1 = linspace(0, sqrt(size_x^2 + size_y^2));
y1 = polyval(p,x1);
plot(x1,y1)
grid minor
set(findall(gcf,'-property','FontSize'),'FontSize',12)
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,strcat(output_path, 'correlation_over_distnace'),'-dpdf','-r0')
hold off


%% Plot filtering results

% Colors for plotting
color_list = {[0 0.4470 0.7410]; [0.8500 0.3250 0.0980]; [0.4660 0.6740 0.1880]; [0.6350 0.0780 0.1840]};

x_axis_scale = 1:new_amount_filtering;
% x_axis_scale = x_axis_scale.^2;

h = figure();
hold on
title(city_name);
xlabel("Emission Resolution [km x km]",'FontWeight','bold');
ylabel("Rel. l_2 error",'FontWeight','bold');
for s=1:amount_SNR
    line = mean(new_error_filter(:,:,s), 2);
    plot(x_axis_scale, line, 'Color', color_list{s},'LineStyle', '--', 'Marker', 'x', 'LineWidth',1.5)
end
for s=1:amount_SNR
    curve1 = prctile(new_error_filter(:,:,s), 80, 2);
    curve2 = prctile(new_error_filter(:,:,s), 20, 2);
    x2 = [x_axis_scale, fliplr(x_axis_scale)];
    inBetween = [curve1', fliplr(curve2')];
    f = fill(x2, inBetween, color_list{s}, 'LineStyle','none');
    set(f,'facealpha',.3)
end
legendStrings = eval(['{' sprintf('''%.0f dB'' ',10*log10(SNR_list)) '}']);
legend(legendStrings);
grid minor
ylim([0 1.1])
xlim([1 x_axis_scale(end)]);
set(findall(gcf,'-property','FontSize'),'FontSize',12)
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,strcat(output_path, 'new_l2_error_over_resolution'),'-dpdf','-r0')
hold off

h = figure();
hold on
title(city_name);
xlabel("Emission Resolution [km x km]",'FontWeight','bold');
ylabel("Rel. l_1 error",'FontWeight','bold');
for s=1:amount_SNR
    line = mean(new_error_filter_l1(:,:,s), 2);
    plot(x_axis_scale, line, 'Color', color_list{s},'LineStyle', '--', 'Marker', 'x', 'LineWidth',1.5)
end
for s=1:amount_SNR
    curve1 = prctile(new_error_filter_l1(:,:,s), 80, 2);
    curve2 = prctile(new_error_filter_l1(:,:,s), 20, 2);
    x2 = [x_axis_scale, fliplr(x_axis_scale)];
    inBetween = [curve1', fliplr(curve2')];
    f = fill(x2, inBetween, color_list{s}, 'LineStyle','none');
    set(f,'facealpha',.3)
end
legendStrings = eval(['{' sprintf('''%.0f dB'' ',10*log10(SNR_list)) '}']);
legend(legendStrings);
ylim([0 1.5])
xlim([1 x_axis_scale(end)]);
grid minor
set(findall(gcf,'-property','FontSize'),'FontSize',12)
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,strcat(output_path, 'new_l1_error_over_resolution'),'-dpdf','-r0')
hold off

%% Plot filtering results DWT

% Colors for plotting
color_list = {[0 0.4470 0.7410]; [0.8500 0.3250 0.0980]; [0.4660 0.6740 0.1880]; [0.6350 0.0780 0.1840]};

x_axis_scale = 1:new_amount_filtering;
% x_axis_scale = x_axis_scale.^2;

h = figure();
hold on
title(city_name);
xlabel("Emission Resolution [km x km]",'FontWeight','bold');
ylabel("Rel. l_2 error",'FontWeight','bold');
for s=1:amount_SNR
    line = mean(new_error_filter_dwt(:,:,s), 2);
    plot(x_axis_scale, line, 'Color', color_list{s},'LineStyle', '--', 'Marker', 'x', 'LineWidth',1.5)
end
for s=1:amount_SNR
    curve1 = prctile(new_error_filter_dwt(:,:,s), 80, 2);
    curve2 = prctile(new_error_filter_dwt(:,:,s), 20, 2);
    x2 = [x_axis_scale, fliplr(x_axis_scale)];
    inBetween = [curve1', fliplr(curve2')];
    f = fill(x2, inBetween, color_list{s}, 'LineStyle','none');
    set(f,'facealpha',.3)
end
legendStrings = eval(['{' sprintf('''%.0f dB'' ',10*log10(SNR_list)) '}']);
legend(legendStrings);
grid minor
ylim([0 1.1])
xlim([1 x_axis_scale(end)]);
set(findall(gcf,'-property','FontSize'),'FontSize',12)
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,strcat(output_path, 'new_l2_error_over_resolution_dwt'),'-dpdf','-r0')
hold off

h = figure();
hold on
title(city_name);
xlabel("Emission Resolution [km x km]",'FontWeight','bold');
ylabel("Rel. l_1 error",'FontWeight','bold');
for s=1:amount_SNR
    line = mean(new_error_filter_l1_dwt(:,:,s), 2);
    plot(x_axis_scale, line, 'Color', color_list{s},'LineStyle', '--', 'Marker', 'x', 'LineWidth',1.5)
end
for s=1:amount_SNR
    curve1 = prctile(new_error_filter_l1_dwt(:,:,s), 80, 2);
    curve2 = prctile(new_error_filter_l1_dwt(:,:,s), 20, 2);
    x2 = [x_axis_scale, fliplr(x_axis_scale)];
    inBetween = [curve1', fliplr(curve2')];
    f = fill(x2, inBetween, color_list{s}, 'LineStyle','none');
    set(f,'facealpha',.3)
end
legendStrings = eval(['{' sprintf('''%.0f dB'' ',10*log10(SNR_list)) '}']);
legend(legendStrings);
grid minor
ylim([0 1.5])
xlim([1 x_axis_scale(end)]);
set(findall(gcf,'-property','FontSize'),'FontSize',12)
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,strcat(output_path, 'new_l1_error_over_resolution_dwt'),'-dpdf','-r0')
hold off

%% Plot filtering results LSQ

% Colors for plotting
color_list = {[0 0.4470 0.7410]; [0.8500 0.3250 0.0980]; [0.4660 0.6740 0.1880]; [0.6350 0.0780 0.1840]};

x_axis_scale = 1:new_amount_filtering;
% x_axis_scale = x_axis_scale.^2;

h = figure();
hold on
title(city_name);
xlabel("Emission Resolution [km x km]",'FontWeight','bold');
ylabel("Rel. l_2 error",'FontWeight','bold');
for s=1:amount_SNR
    line = mean(new_error_filter_lsq(:,:,s), 2);
    plot(x_axis_scale, line, 'Color', color_list{s},'LineStyle', '--', 'Marker', 'x', 'LineWidth',1.5)
end
for s=1:amount_SNR
    curve1 = prctile(new_error_filter_lsq(:,:,s), 80, 2);
    curve2 = prctile(new_error_filter_lsq(:,:,s), 20, 2);
    x2 = [x_axis_scale, fliplr(x_axis_scale)];
    inBetween = [curve1', fliplr(curve2')];
    f = fill(x2, inBetween, color_list{s}, 'LineStyle','none');
    set(f,'facealpha',.3)
end
legendStrings = eval(['{' sprintf('''%.0f dB'' ',10*log10(SNR_list)) '}']);
legend(legendStrings);
grid minor
ylim([0 1.1])
xlim([1 x_axis_scale(end)]);
set(findall(gcf,'-property','FontSize'),'FontSize',12)
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,strcat(output_path, 'new_l2_error_over_resolution_lsq'),'-dpdf','-r0')
hold off

h = figure();
hold on
title(city_name);
xlabel("Emission Resolution [km x km]",'FontWeight','bold');
ylabel("Rel. l_1 error",'FontWeight','bold');
for s=1:amount_SNR
    line = mean(new_error_filter_l1_lsq(:,:,s), 2);
    plot(x_axis_scale, line, 'Color', color_list{s},'LineStyle', '--', 'Marker', 'x', 'LineWidth',1.5)
end
for s=1:amount_SNR
    curve1 = prctile(new_error_filter_l1_lsq(:,:,s), 80, 2);
    curve2 = prctile(new_error_filter_l1_lsq(:,:,s), 20, 2);
    x2 = [x_axis_scale, fliplr(x_axis_scale)];
    inBetween = [curve1', fliplr(curve2')];
    f = fill(x2, inBetween, color_list{s}, 'LineStyle','none');
    set(f,'facealpha',.3)
end
legendStrings = eval(['{' sprintf('''%.0f dB'' ',10*log10(SNR_list)) '}']);
legend(legendStrings);
grid minor
ylim([0 1.5])
xlim([1 x_axis_scale(end)]);
set(findall(gcf,'-property','FontSize'),'FontSize',12)
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,strcat(output_path, 'new_l1_error_over_resolution_lsq'),'-dpdf','-r0')
hold off


%% Plot filtering results combined
% Colors for plotting
color_list = {[0, 0.4470, 0.7410]; [0.8500, 0.3250, 0.0980]; [0.9290, 0.6940, 0.1250]; [0.6350 0.0780 0.1840]};

x_axis_scale = 1:new_amount_filtering;
% x_axis_scale = x_axis_scale.^2;

h = figure();
hold on
title(city_name);
xlabel("Emission Resolution [km x km]",'FontWeight','bold');
ylabel("Rel. l_2 error",'FontWeight','bold');
for s=1:amount_SNR
    line = mean(new_error_filter(:,:,s), 2);
    plot(x_axis_scale, line, 'Color', color_list{1},'LineStyle', '--', 'Marker', 'x', 'LineWidth',1.5)
    line = mean(new_error_filter_dwt(:,:,s), 2);
    plot(x_axis_scale, line, 'Color', color_list{2},'LineStyle', '--', 'Marker', 'x', 'LineWidth',1.5)
    line = mean(new_error_filter_lsq(:,:,s), 2);
    plot(x_axis_scale, line, 'Color', color_list{3},'LineStyle', '--', 'Marker', 'x', 'LineWidth',1.5)
end
% for s=1:amount_SNR
%     curve1 = prctile(new_error_filter(:,:,s), 80, 2);
%     curve2 = prctile(new_error_filter(:,:,s), 20, 2);
%     x2 = [x_axis_scale, fliplr(x_axis_scale)];
%     inBetween = [curve1', fliplr(curve2')];
%     f = fill(x2, inBetween, color_list{1}, 'LineStyle','none');
%     set(f,'facealpha',.3)
%     curve1 = prctile(new_error_filter_dwt(:,:,s), 80, 2);
%     curve2 = prctile(new_error_filter_dwt(:,:,s), 20, 2);
%     x2 = [x_axis_scale, fliplr(x_axis_scale)];
%     inBetween = [curve1', fliplr(curve2')];
%     f = fill(x2, inBetween, color_list{2}, 'LineStyle','none');
%     set(f,'facealpha',.3)
%     curve1 = prctile(new_error_filter_lsq(:,:,s), 80, 2);
%     curve2 = prctile(new_error_filter_lsq(:,:,s), 20, 2);
%     x2 = [x_axis_scale, fliplr(x_axis_scale)];
%     inBetween = [curve1', fliplr(curve2')];
%     f = fill(x2, inBetween, color_list{3}, 'LineStyle','none');
%     set(f,'facealpha',.3)
% end
% legendStrings = eval(['{' sprintf('''%.0f dB'' ',10*log10(SNR_list)) '}']);
legend("SR", "SR DWT", "LS");
grid minor
ylim([0 1.0])
xlim([1 x_axis_scale(end)]);
set(findall(gcf,'-property','FontSize'),'FontSize',12)
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,strcat(output_path, 'new_l2_error_over_resolution_combined'),'-dpdf','-r0')
hold off

%% Plot reconstruction pathway
x_axis_labels = {};
for k=1:amount_captures
    j = floor(captured_points(k));
    x_axis_labels{end + 1} = sprintf('%d - %d %%',0,floor((j/(size_x * size_y))*100));
end
x_axis_labels = categorical(x_axis_labels, x_axis_labels, 'Ordinal',true);


h = figure();
hold on
bar(x_axis_labels, [captured_mean(:,1) captured_mean(:,2) captured_mean(:,3)], 'linewidth',1.5);
xlabel("Highest emissions inventory",'FontWeight','bold');
ylabel("Percent captured by reconstruction",'FontWeight','bold');
ylim([0 1]);
lgd = legend(legendStrings);
lgd.Location = 'northwest';
grid minor
%set(gca, 'YScale', 'log')
set(findall(gcf,'-property','FontSize'),'FontSize',12)
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,strcat(output_path, 'emission_reconstruction_pathway_percent'),'-dpdf','-r0')

%% Plot reconstruction pathway DWT

h = figure();
hold on
bar(x_axis_labels, [captured_dwt_mean(:,1) captured_dwt_mean(:,2) captured_dwt_mean(:,3)], 'linewidth',1.5);
xlabel("Highest emissions inventory",'FontWeight','bold');
ylabel("Percent captured by reconstruction",'FontWeight','bold');
ylim([0 1]);
lgd = legend(legendStrings);
lgd.Location = 'northwest';
grid minor
%set(gca, 'YScale', 'log')
set(findall(gcf,'-property','FontSize'),'FontSize',12)
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,strcat(output_path, 'emission_dwt_reconstruction_pathway_percent'),'-dpdf','-r0')