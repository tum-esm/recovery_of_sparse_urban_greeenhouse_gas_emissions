%% Description

% Anaylzes a sensing matrix towards its compressed sensing capabilities,
% using the incoherence parameter and by performing a Monte Carlo
% simulation to determine a lower bound on the RIP.

%% Load gurobi as cvx solver
cvx_solver Gurobi

%% functions
get2 = @(x, i) reshape(x(i, :), size(x, 2), 1);
get3 = @(x, i) reshape(x(i, :, :), size(x, 2), size(x, 3));
get2of3 = @(x, i, j) reshape(x(i, j, :), size(x, 3));


%% Load sensing matrices
sensing_matrix_path = "../data/footprint/synthtetic.nc";


%% set variables
sensing_threshhold = 1e-9;
dwt_level = 3;

% domain size for footprints
size_x = 32;
size_y = 32;

% name normal stations
measurement_station_names = ["A"; "B"; "C"; "D"; "E"; "F"; "G"; "H"; "I"; "J";"K";"L";"M";"N";"O";"P"];

% name highest wind coverage in wind variation footprint
% measurement_station_names = ["A_6.283185307179585"; "B_6.283185307179585"; "C_6.283185307179585"; "D_6.283185307179585"; 
%                         "E_6.283185307179585"; "F_6.283185307179585"; "G_6.283185307179585"; "H_6.283185307179585"; "I_6.283185307179585"; 
%                         "J_6.283185307179585";"K_6.283185307179585";"L_6.283185307179585";"M_6.283185307179585";"N_6.283185307179585";
%                         "O_6.283185307179585";"P_6.283185307179585"];

% name lowest wind coverage in wind variation footprint
% measurement_station_names = ["A_0.41887902047863906"; "B_0.41887902047863906"; "C_0.41887902047863906"; "D_0.41887902047863906"; 
%                         "E_0.41887902047863906"; "F_0.41887902047863906"; "G_0.41887902047863906"; "H_0.41887902047863906"; "I_0.41887902047863906"; 
%                         "J_0.41887902047863906";"K_0.41887902047863906";"L_0.41887902047863906";"M_0.41887902047863906";"N_0.41887902047863906";
%                         "O_0.41887902047863906";"P_0.41887902047863906"];

amount_measurement_stations = size(measurement_station_names, 1);

percent_measurement_stations_to_use = 0.015;

amount_measurement_stations_to_use = floor(size_x * size_y * percent_measurement_stations_to_use);

if amount_measurement_stations_to_use > amount_measurement_stations
    amount_measurement_stations_to_use  = amount_measurement_stations;
end


%% Prepare Sensing matrix

collection_sensing_matrix = load_collection_footprint(sensing_matrix_path, measurement_station_names, size_x, size_y);

% create the sensing matrix
sensing_matrix = [];
for j=1:amount_measurement_stations_to_use
    tmp = get3(collection_sensing_matrix, j);
    sensing_matrix = [sensing_matrix; tmp];
end
amount_measurements = size(sensing_matrix, 1);

%% Apply threshholding to sensing matrix
temp_sensing_mat = sensing_matrix;
temp_sensing_mat(temp_sensing_mat < sensing_threshhold) = 0;
reduction_mask = ~all(~temp_sensing_mat,1);

T = eye(size_x * size_y);
T = T(reduction_mask, :);

sensing_matrix_tilde = sensing_matrix * transpose(T);

%% Create Wavelet sensing matrix

[dwt_sensing_mat, S, wave_size] = DWT_matrix_rows(sensing_matrix, size_x, size_y, dwt_level);

T_dwt = eye(wave_size);

sensing_matrix_tilde_dwt = dwt_sensing_mat * transpose(T_dwt);


%% Calculate the coherence criteria
remove_diag = @(C)C-diag(diag(C));
Cor = @(Phi)remove_diag(abs(Phi'*Phi));
normalize_columns = @(X)X/(diag(sqrt(diag(X'*X))));
maxall = @(C)max(C(:));
meanall = @(C)mean(C(:));
mu = @(Phi)maxall(Cor(normalize_columns(Phi)));

fprintf("coherence of sensing matrix is %f\n", mu(sensing_matrix_tilde));
fprintf("coherence of sensing matrix in the wavelet domain is %f\n", mu(sensing_matrix_tilde_dwt));

%% Monte Carlo for Restricted Isometry Property

Phi = normalize_columns(sensing_matrix_tilde);

N = size(Phi, 2); % amount DoF
m = size(Phi, 1); % amount measurmenets

s = 15;
amount_iterations = 2000;
t = linspace(0,5,100);
v = [];
for i=1:amount_iterations
    I = randperm(N); I = I(1:s);
    v = [v; svd(Phi(:, I))];
end
figure();
h = hist(v,t);
bar(t,h);
% robust lines
xline(1-sqrt(2)+1, 'r');
xline(1+sqrt(2)-1, 'r');
% noiseless lines
xline(0, 'r');
xline(2, 'r');
title(sprintf('n=%d, m=%d, s=%d',N, m,s ));

%% Monte Carlo for Restricted Isometry Property DWT

Phi = normalize_columns(sensing_matrix_tilde_dwt);

N = size(Phi, 2); % amount DoF
m = size(Phi, 1); % amount measurmenets

s = 15;
amount_iterations = 2000;
t = linspace(0,5,100);
v = [];
for i=1:amount_iterations
    I = randperm(N); I = I(1:s);
    v = [v; svd(Phi(:, I))];
end
figure();
h = hist(v,t);
bar(t,h);
% robust lines
xline(1-sqrt(2)+1, 'r');
xline(1+sqrt(2)-1, 'r');
% noiseless lines
xline(0, 'r');
xline(2, 'r');
title(sprintf('n=%d, m=%d, s=%d',N, m,s ));