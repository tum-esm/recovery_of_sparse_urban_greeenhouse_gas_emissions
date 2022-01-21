%% Description

% Does the same as "city_recovery.m", with the difference that noise is
% considered.

%% Load gurobi as cvx solver
cvx_solver Gurobi

%% functions
get2 = @(x, i) reshape(x(i, :), size(x, 2), 1);
get3 = @(x, i) reshape(x(i, :, :), size(x, 2), size(x, 3));
get2of3 = @(x, i, j) reshape(x(i, j, :), size(x, 3));
normalize_columns = @(X)X/(diag(sqrt(diag(X'*X))));
get_normalization_vector = @(X)eye(size(X, 2))/diag(sqrt(diag(X'*X)));
logn = @(a, n) log(a) / log(n);


%% Input data

map = "co2_vienna.nc";
city_name = "Vienna";

SNR = 100;

output_path_folder = '../output/case_studies';

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
percent_measurement_stations_to_use = 0.03;

ncfile = "../data/emission_map/" + map;
co2 = ncread(ncfile, 'CO2Diffuse');
co2 = 1e-4 * co2; % convert co2 scale

%% Preprocessing

size_x = size(co2, 1);
size_y = size(co2, 2);
co2_vec = reshape(co2, size_x * size_y, 1);




collection_sensing_matrix = load_collection_footprint(sensing_matrix_path, measurement_station_names, size_x, size_y);


amount_measurement_stations_to_use = floor(size_x * size_y * percent_measurement_stations_to_use);
if amount_measurement_stations_to_use > amount_measurement_stations
    error('not enough measurement stations in the footprint!');
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


%% Create Emission Map
emission_map = reshape(co2_vec, size_x, size_y);
emission_map_vec = reshape(emission_map, size_x * size_y, 1);

real_s = sum(emission_map_vec ~= 0);


%% Create observation
obs_vec = double(sensing_matrix * emission_map_vec);

%% Create noise
amount_observations = size(obs_vec, 1);
signal_power = sum_square(obs_vec) / amount_observations;
sigma_noise = (signal_power/SNR)^0.5;
e = sigma_noise * randn(amount_observations, 1);

obs_vec_noiseless = obs_vec;
obs_vec = obs_vec + e;

%% CS reconstruction
emission_L1_vec = optimizeL1_noise(sensing_matrix_tilde, obs_vec, sum_square(e));
emission_L1_vec = transpose(T) * emission_L1_vec;


emission_L1 = reshape(emission_L1_vec, size_x, size_y);

%% Wavelet CS
emission_L1_dwt_vec = optimizeL1_noise(sensing_matrix_tilde_dwt, obs_vec, sum_square(e));

emission_L1_dwt_vec = transpose(T_dwt) * emission_L1_dwt_vec;

emission_L1_by_DWT = waverec2(emission_L1_dwt_vec, S, 'haar');

emission_L1_vec_by_DWT = reshape(emission_L1_by_DWT, size_x * size_y, 1);
clear emission_L1_dwt_vec

%% Least Square
emission_l2_vec = optimizeL2_noise(sensing_matrix_tilde, obs_vec, sum_square(e));
emission_l2_vec = transpose(T) * emission_l2_vec;

emission_l2 = reshape(emission_l2_vec, size_x, size_y);


%% Generate plotting path
output_path = strcat(output_path_folder, "/", city_name, "/", sprintf("%.1fdB_",10*log10(SNR)),"noise_", "s_", string(real_s), "ms_", string(percent_measurement_stations_to_use), "/");

mkdir(output_path);

%% Printing information
fileID = fopen(strcat(output_path, 'info.txt'),'w');
fprintf(fileID, "Information file for %s\n", city_name);
fprintf(fileID, "\n");
fprintf(fileID, "n is %d\n", size_x * size_y);
fprintf(fileID, "m is %d\n", size(sensing_matrix_tilde, 1));
fprintf(fileID, "Amount measurement stations used is %d\n", amount_measurement_stations_to_use);
fprintf(fileID, "SNR is %.1fdB\n", 10*log10(SNR));
fprintf(fileID, "\n");
fprintf(fileID, "The l2 reconstruction error using sparse reconstruction is %f\n", norm(emission_L1_vec - emission_map_vec, 2)/norm(emission_map_vec,2));
fprintf(fileID, "The l1 reconstruction error using sparse reconstruction is %f\n", norm(emission_L1_vec - emission_map_vec, 1)/norm(emission_map_vec,1));
fprintf(fileID, "\n");
fprintf(fileID, "The l2 reconstruction error using sparse DWT reconstruction is %f\n", norm(emission_L1_vec_by_DWT - emission_map_vec, 2)/norm(emission_map_vec,2));
fprintf(fileID, "The l1 reconstruction error using sparse DWT reconstruction is %f\n", norm(emission_L1_vec_by_DWT - emission_map_vec, 1)/norm(emission_map_vec,1));
fprintf(fileID, "\n");
fprintf(fileID, "The l2 reconstruction error using least squares is %f\n", norm(emission_l2_vec - emission_map_vec, 2)/norm(emission_map_vec,2));
fprintf(fileID, "The l1 reconstruction error using least squares is %f\n", norm(emission_l2_vec - emission_map_vec, 1)/norm(emission_map_vec,1));
fprintf(fileID, "\n");
%% Plotting l1 reconstructed emission map vs real emission map
max_color = max([max(max(emission_map)), max(max(emission_L1)), max(max(emission_L1_by_DWT))]);
h=figure('Position', [100, 100, 3*512, 512]);
ax = gca;
h1 = subplot(1,3,1);
imagesc(emission_L1, [0, max_color]);
ax1 = gca;
title("reconstructed emission");
hold on
h2 = subplot(1,3,2);
imagesc(emission_L1_by_DWT, [0, max_color]);
ax2 = gca;
title("reconstructed emission DWT");
h3 = subplot(1,3,3);
imagesc(emission_map, [0, max_color]);
ax3 = gca;
title("real emission");
sgtitle("Comparison L1 reconstruction to Real Emissions")
%grid minor
set(findall(gcf,'-property','FontSize'),'FontSize',12)
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,strcat(output_path, 'maps_l1_real_comparision'),'-dpdf','-r0')


%% Plotting elemenwtise l2 error in reconstrucion
h=figure();
imagesc(sqrt((emission_L1 - emission_map).*(emission_L1 - emission_map))/norm(emission_map_vec, 2));
hold on
grid minor
colorbar
set(findall(gcf,'-property','FontSize'),'FontSize',12)
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
title("Rel. L2 error for SR");
print(h,strcat(output_path, 'elementwise_l2_error_BP'),'-dpdf','-r0')

%% Plotting elemenwtise l2 error in reconstrucion for DWT
h=figure();
imagesc(sqrt((emission_L1_by_DWT - emission_map).*(emission_L1_by_DWT - emission_map))/norm(emission_map_vec, 2));
hold on
grid minor
colorbar
set(findall(gcf,'-property','FontSize'),'FontSize',12)
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
title("Rel. L2 error for SR with DWT");
print(h,strcat(output_path, 'elementwise_l2_error_DWT'),'-dpdf','-r0')

%% Plotting elemenwtise l2 error in reconstrucion for least squares
h=figure();
imagesc(sqrt((emission_l2 - emission_map).*(emission_l2 - emission_map))/norm(emission_map_vec, 2));
hold on
grid minor
colorbar
set(findall(gcf,'-property','FontSize'),'FontSize',12)
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
title("Rel. L2 error for LS");
print(h,strcat(output_path, 'elementwise_l2_error_least_square'),'-dpdf','-r0')


%% Plotting Smoothed errors
max_smoothing = max(size_x, size_y);

l1_smoothing_l2_error = zeros(max_smoothing,1);
l1dwt_smoothing_l2_error = zeros(max_smoothing,1);
l2_smoothing_l2_error = zeros(max_smoothing,1);

l1_smoothing_l2_error(1) = norm(emission_L1_vec - emission_map_vec, 2)/norm(emission_map_vec,2);
l1dwt_smoothing_l2_error(1) = norm(emission_L1_vec_by_DWT - emission_map_vec, 2)/norm(emission_map_vec,2);
l2_smoothing_l2_error(1) = norm(emission_l2_vec - emission_map_vec, 2)/norm(emission_map_vec,2);
for j=2:max_smoothing
    smoothed_l1 = smoothing(emission_L1, j);
    smoothed_l1_dwt = smoothing(emission_L1_by_DWT, j);
    smoothed_l2 = smoothing(emission_l2, j);
    smoothed_real = smoothing(emission_map, j);
    
    smoothed_l1 = reshape(smoothed_l1, size_x * size_y, 1);
    smoothed_l1_dwt = reshape(smoothed_l1_dwt, size_x * size_y, 1);
    smoothed_l2 = reshape(smoothed_l2, size_x * size_y, 1);
    smoothed_real = reshape(smoothed_real, size_x * size_y, 1);
    
    l1_smoothing_l2_error(j) = norm(smoothed_l1 - smoothed_real, 2)/norm(smoothed_real,2);
    l1dwt_smoothing_l2_error(j) = norm(smoothed_l1_dwt - smoothed_real, 2)/norm(smoothed_real,2);
    l2_smoothing_l2_error(j) = norm(smoothed_l2 - smoothed_real, 2)/norm(smoothed_real,2);
end

x_axis_scale = 1:max_smoothing;

h = figure();
hold on
xlabel("Emission Resolution [km x km]",'FontWeight','bold');
ylabel("Rel. l_2 error",'FontWeight','bold');
plot(x_axis_scale, l1_smoothing_l2_error,'LineStyle', '--', 'Marker', 'x', 'LineWidth',1.5)
plot(x_axis_scale, l1dwt_smoothing_l2_error,'LineStyle', '--', 'Marker', 'x', 'LineWidth',1.5)
plot(x_axis_scale, l2_smoothing_l2_error,'LineStyle', '--', 'Marker', 'x', 'LineWidth',1.5)
legend("SR", "SR DWT", "LS");
grid minor
ylim([0 inf])
xlim([1 x_axis_scale(end)]);
set(findall(gcf,'-property','FontSize'),'FontSize',12)
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,strcat(output_path, 'l2_error_over_resolution'),'-dpdf','-r0')
hold off

%% Write 5 km x 5 km smoothing error to file
fprintf(fileID, "l2 5 km x 5 km smoothed error for CS is %f\n", l1_smoothing_l2_error(5));
fprintf(fileID, "l2 5 km x 5 km smoothed error for CS DWT is %f\n", l1dwt_smoothing_l2_error(5));
fprintf(fileID, "l2 5 km x 5 km smoothed error for LSQ is %f\n", l2_smoothing_l2_error(5));
fprintf(fileID, "\n");
fprintf(fileID, "l1 5 km x 5 km smoothed error for CS is %f\n", l1_smoothing_l1_error(5));
fprintf(fileID, "l1 5 km x 5 km smoothed error for CS DWT is %f\n", l1dwt_smoothing_l1_error(5));
fprintf(fileID, "l1 5 km x 5 km smoothed error for LSQ is %f\n", l2_smoothing_l1_error(5));
fprintf(fileID, "\n");
%% Write total error to file
total_error_l1 = (sum(emission_L1_vec) - sum(emission_map_vec))/sum(emission_map_vec);
total_error_l1_dwt = (sum(emission_L1_vec_by_DWT) - sum(emission_map_vec))/sum(emission_map_vec);
total_error_lsq = (sum(emission_l2_vec) - sum(emission_map_vec))/sum(emission_map_vec);

fprintf(fileID, "Total emission error using sparse reconstruction is %f\n", total_error_l1);
fprintf(fileID, "Total emission error using sparse reconstruction with DWT is %f\n", total_error_l1_dwt);
fprintf(fileID, "Total emission error using least square method is %f\n", total_error_lsq);
fprintf(fileID, "\n");

%% Plot Emission strengths

[map_sorted, idx_map] = sort(emission_map_vec, 'descend');
[l1_sorted, idx_l1] = sort(emission_L1_vec, 'descend');
[l1_dwt_sorted, idx_l1_dwt] = sort(emission_L1_vec_by_DWT, 'descend');
[l2_sorted, idx_l2] = sort(emission_l2_vec, 'descend');

h = figure();
plot(map_sorted, 'linewidth',1.5);
%title("Distribution of Emissions");
hold on
plot(l1_sorted, 'linewidth',1.5);
plot(l1_dwt_sorted, 'linewidth',1.5);
plot(l2_sorted, 'linewidth',1.5);
%xlabel("Rel. sparsity s / DoF",'FontWeight','bold');
ylabel("Emission strength",'FontWeight','bold');
lgd = legend('Inventory', 'SR', 'SR DWT','LS');
lgd.Location = 'southwest';
grid minor
set(gca, 'YScale', 'log')
set(findall(gcf,'-property','FontSize'),'FontSize',12)
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,strcat(output_path, 'distribution_recovered_emissions'),'-dpdf','-r0')

%% Give Amount of Same highest Emissions
perfect_path = [];
CS_path = [];
CS_DWT_path = [];
LS_path = [];
for j=1:(size_x * size_y)
    equal_same = length(intersect(idx_map(1:j), idx_map(1:j)));
    equal_l1 = length(intersect(idx_map(1:j), idx_l1(1:j)));
    equal_l1_dwt = length(intersect(idx_map(1:j), idx_l1_dwt(1:j)));
    equal_l2 = length(intersect(idx_map(1:j), idx_l2(1:j)));
    
    perfect_path = [perfect_path; equal_same];
    CS_path = [CS_path; equal_l1];
    CS_DWT_path = [CS_DWT_path; equal_l1_dwt];
    LS_path = [LS_path; equal_l2];
end


h = figure();
plot(perfect_path, 'linewidth',1.5);
%title("Reconstruction Path");
hold on
plot(CS_path, 'linewidth',1.5);
plot(CS_DWT_path, 'linewidth',1.5);
plot(LS_path, 'linewidth',1.5);
xlabel("n^{th} Highest emissions inventory",'FontWeight','bold');
ylabel("n^{th} Highest emissions reconstruction",'FontWeight','bold');
lgd = legend('Ideal', 'SR', 'SR DWT','LS');
lgd.Location = 'northwest';
grid minor
% set(gca, 'YScale', 'log')
set(findall(gcf,'-property','FontSize'),'FontSize',12)
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,strcat(output_path, 'emission_reconstruction_pathway'),'-dpdf','-r0')

%% Reconstruction pathway percent
perfect_path = [];
CS_path = [];
CS_DWT_path = [];
LS_path = [];
x_axis_labels = {};
x_axis = linspace(nthroot(0, 2.5), nthroot(1, 2.5), 10).^2.5;
for k=1:8
    j = floor(x_axis(k) * size_x * size_y);
    if j > 1
        equal_same = length(intersect(idx_map(1:j), idx_map(1:j)));
        equal_l1 = length(intersect(idx_map(1:j), idx_l1(1:j))) / equal_same;
        equal_l1_dwt = length(intersect(idx_map(1:j), idx_l1_dwt(1:j))) / equal_same;
        equal_l2 = length(intersect(idx_map(1:j), idx_l2(1:j))) / equal_same;

        x_axis_labels{end + 1} = sprintf('%d - %.1f %%',0,(j/real_s)*100);
        perfect_path = [perfect_path; equal_same];
        CS_path = [CS_path; equal_l1];
        CS_DWT_path = [CS_DWT_path; equal_l1_dwt];
        LS_path = [LS_path; equal_l2];
    end
end
x_axis_labels = categorical(x_axis_labels, x_axis_labels, 'Ordinal',true);

h = figure();
title("Qualitative Measure");
hold on
bar(x_axis_labels, [CS_path CS_DWT_path LS_path], 'linewidth',1.5);
xlabel("Highest emissions inventory",'FontWeight','bold');
ylabel("Fraction reconstructed",'FontWeight','bold');
lgd = legend('SR', 'SR DWT','LS');
lgd.Location = 'southeast';
ylim([0 1])
grid minor
%set(gca, 'YScale', 'log')
set(findall(gcf,'-property','FontSize'),'FontSize',12)
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,strcat(output_path, 'emission_reconstruction_pathway_percent'),'-dpdf','-r0')

%% Error path of highest emissions
CS_path = [];
CS_DWT_path = [];
LS_path = [];
for j=1:(size_x * size_y)
    equal_l1 = norm(emission_map_vec(idx_map(1:j)) - emission_L1_vec(idx_map(1:j)),2)/norm(emission_map_vec(idx_map(1:j)),2);
    equal_l1_dwt = norm(emission_map_vec(idx_map(1:j)) - emission_L1_vec_by_DWT(idx_map(1:j)),2)/norm(emission_map_vec(idx_map(1:j)),2);
    equal_l2 = norm(emission_map_vec(idx_map(1:j)) - emission_l2_vec(idx_map(1:j)),2)/norm(emission_map_vec(idx_map(1:j)),2);

    CS_path = [CS_path; equal_l1];
    CS_DWT_path = [CS_DWT_path; equal_l1_dwt];
    LS_path = [LS_path; equal_l2];
end


h = figure();
plot(CS_path, 'linewidth',1.5);
%title("Reconstruction Error Path");
hold on
plot(CS_DWT_path, 'linewidth',1.5);
plot(LS_path, 'linewidth',1.5);
xlabel("n^{th} Highest emissions inventory",'FontWeight','bold');
ylabel("Rel. l_2 error",'FontWeight','bold');
lgd = legend('SR', 'SR DWT','LS');
lgd.Location = 'northwest';
grid minor
% set(gca, 'YScale', 'log')
set(findall(gcf,'-property','FontSize'),'FontSize',12)
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,strcat(output_path, 'reconstruction_error_pathway'),'-dpdf','-r0')

%% Average of relative errors of clustered emissions
CS_path = [];
CS_DWT_path = [];
LS_path = [];
x_axis_labels = {};
j2 = 1;
for k=1:((size_x * size_y)-1)/10:(size_x * size_y)
    j = floor(k);
    if j ~= 1
        equal_l1 = mean(abs(emission_map_vec(idx_map(j2:j)) - emission_L1_vec(idx_map(j2:j))) ./ emission_map_vec(idx_map(j2:j)));
        equal_l1_dwt = mean(abs(emission_map_vec(idx_map(j2:j)) - emission_L1_vec_by_DWT(idx_map(j2:j))) ./ emission_map_vec(idx_map(j2:j)));
        equal_l2 = mean(abs(emission_map_vec(idx_map(j2:j)) - emission_l2_vec(idx_map(j2:j))) ./ emission_map_vec(idx_map(j2:j)));
        
        x_axis_labels{end + 1} = sprintf('%d - %d %%',floor((j2/(size_x * size_y))*100),floor((j/(size_x * size_y))*100));
        CS_path = [CS_path; equal_l1];
        CS_DWT_path = [CS_DWT_path; equal_l1_dwt];
        LS_path = [LS_path; equal_l2];
        
        j2 = j;
    end
end

x_axis_labels = categorical(x_axis_labels, x_axis_labels, 'Ordinal',true);
h = figure();
bar(x_axis_labels, [CS_path CS_DWT_path LS_path], 'linewidth',1.5);
%title("Reconstruction Error Path");
hold on
%bar(x_axis_labels, CS_DWT_path, 'linewidth',1.5);
%bar(x_axis_labels, LS_path, 'linewidth',1.5);
xlabel("Highest emissions in inventory",'FontWeight','bold');
ylabel("Mean rel. error",'FontWeight','bold');
lgd = legend('SR', 'SR DWT','LS');
lgd.Location = 'northwest';
grid minor
set(gca, 'YScale', 'log')
set(findall(gcf,'-property','FontSize'),'FontSize',12)
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,strcat(output_path, 'reconstruction_clustered_mean_error'),'-dpdf','-r0')

%% Average of relative errors of clustered emissions, logarithmic
CS_path = [];
CS_DWT_path = [];
LS_path = [];
x_axis_labels = {};
j2 = 1;
for k=1:8
    j = floor(x_axis(k) * size_x * size_y);
    if j > 1
        equal_l1 = mean(abs(emission_map_vec(idx_map(j2:j)) - emission_L1_vec(idx_map(j2:j))) ./ emission_map_vec(idx_map(j2:j)));
        equal_l1_dwt = mean(abs(emission_map_vec(idx_map(j2:j)) - emission_L1_vec_by_DWT(idx_map(j2:j))) ./ emission_map_vec(idx_map(j2:j)));
        equal_l2 = mean(abs(emission_map_vec(idx_map(j2:j)) - emission_l2_vec(idx_map(j2:j))) ./ emission_map_vec(idx_map(j2:j)));
        
        if j2 == 1
            j2 = 0;
        end
        x_axis_labels{end + 1} = sprintf('%.1f - %.1f %%',(j2/(size_x * size_y))*100,(j/(size_x * size_y))*100);
        CS_path = [CS_path; equal_l1];
        CS_DWT_path = [CS_DWT_path; equal_l1_dwt];
        LS_path = [LS_path; equal_l2];
        
        j2 = j;
    end
end

x_axis_labels = categorical(x_axis_labels, x_axis_labels, 'Ordinal',true);
h = figure();
bar(x_axis_labels, [CS_path CS_DWT_path LS_path], 'linewidth',1.5);
ylim([0.006 1.9])
title("Quantitative Measure");
hold on
%bar(x_axis_labels, CS_DWT_path, 'linewidth',1.5);
%bar(x_axis_labels, LS_path, 'linewidth',1.5);
xlabel("Highest emissions in inventory",'FontWeight','bold');
ylabel("Mean rel. error",'FontWeight','bold');
lgd = legend('SR', 'SR DWT','LS');
lgd.Location = 'southeast';
grid minor
set(gca, 'YScale', 'log')
set(findall(gcf,'-property','FontSize'),'FontSize',12)
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,strcat(output_path, 'reconstruction_logarithmic_clustered_mean_error'),'-dpdf','-r0')

%% Give Amount of Same highest Emissions
equal_l1 = find((idx_map - idx_l1)~=0, 1, 'first') - 1;
equal_l1_dwt = find((idx_map - idx_l1_dwt)~=0, 1, 'first') - 1;
equal_l2 = find((idx_map - idx_l2)~=0, 1, 'first') - 1;

fprintf(fileID, "Emission strength order equal for first %d indices for CS\n", equal_l1);
fprintf(fileID, "Emission strength order equal for first %d indices for CS DWT\n", equal_l1_dwt);
fprintf(fileID, "Emission strength order equal for first %d indices for least square\n", equal_l2);
fprintf(fileID, "\n");
fclose(fileID);

    