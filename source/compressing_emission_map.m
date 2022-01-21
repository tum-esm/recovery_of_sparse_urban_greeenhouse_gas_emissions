%% Description

% Generates a table and xls file which shows the compressability of different
% inventories, both in the spatial and wavelet domain. As compressability
% measures, the Gini index and \sigma_10% is used. (For details, see the
% paper of this repository).

%% Code
% give all the names for the emission maps you want to
% analyze, which have to be located in ../data/emission_map
maps = ["co2_berlin.nc"; 
        "co2_hamburg.nc"; 
        "co2_london.nc";
        "co2_munich.nc";
        "co2_paris.nc";
        "co2_vienna.nc"];

dwt_level = 3;
gini_indices = [];
gini_indices_wavelet = [];

relative_errors = [];
relative_errors_wavelet = [];
L1_relative_errors = [];
L1_relative_errors_wavelet = [];


for i=1:length(maps)
    map = maps(i);
    ncfile = "../data/emission_map/" + map;
    co2 = ncread(ncfile, 'CO2Diffuse');
    
    
    % co2 = imresize(co2, [40, 40]);
    size_x = size(co2, 1);
    size_y = size(co2, 2);
    
    co2 = reshape(co2, size_x * size_y, 1);
    
    n = ceil(size_x * size_y * 0.1);
    
    
    gini_indices = [gini_indices; get_GiniIndex(co2)];
    
    reduced_co2 = get_LargestValueAbsolute(co2, n);
    
    rel_error = norm(co2 - reduced_co2, 2)/norm(co2, 2);
    relative_errors = [relative_errors; rel_error];
    
    
    [co2_dwt, S] = wavedec2(co2, dwt_level, 'haar');
    
    gini_indices_wavelet = [gini_indices_wavelet; get_GiniIndex(co2_dwt)];
    
    reduced_co2_dwt = get_LargestValueAbsolute(co2_dwt, n);
    
    rec_co2 = waverec2(co2_dwt, S, 'haar');
    rec_co2_reduced = waverec2(reduced_co2_dwt, S, 'haar');
    
    rel_error_dwt = norm(co2 - rec_co2_reduced, 2)/norm(co2, 2);
    relative_errors_wavelet = [relative_errors_wavelet; rel_error_dwt];
end

T = table(maps, gini_indices, gini_indices_wavelet, relative_errors, relative_errors_wavelet)
writetable(T, "../output/compressability_table.xls");