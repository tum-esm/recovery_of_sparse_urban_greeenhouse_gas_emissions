function [emission_map] = create_s_sparse_emission_map(s,size_x,size_y,centered)
% create_s_sparse_emission_map  Create a s-sparse emission map.
%   C = create_s_sparse_emission_map(s, size_x, size_y) sparsity s in percent or integer and size x and y of emission map.
%
%   C = create_s_sparse_emission_map(s, size_x, size_y) centered true if
%   sparse emissions should be centered
if s>1
    s = s/(size_x * size_y);
end
if nargin < 4
    centered = false;
end
emission_map = sprandn(size_x, size_y, s);
mask = emission_map~=0;
emission_map(mask) = emission_map(mask) + 10;
emission_map = full(emission_map);

if centered
   emission_map = zeros(size_x, size_y);
    s_start = 0;
    s_needed = s * size_x * size_y;
    while s_start < s_needed

        sigma_x = 0.15 * size_x;
        sigma_y = 0.15 * size_y;

        % probability distribution of emission location
        emission_x = sigma_x * randn() + (size_x/2);
        emission_y = sigma_y * randn() + (size_y/2);
        
        c = max(min(floor(emission_x), size_x), 1);
        d = max(min(floor(emission_y), size_y), 1);

        if emission_map(c,d) == 0
            emission_map(c,d) = 10 + randn();
            s_start = s_start +1;
        end
    end
end

