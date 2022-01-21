function index = get_GiniIndex(x)
    % convert the input matrix into vector if the input x is a matrix
    % Using the formula given in
    % https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=5934357
    if size(x,2)>1
        c = x(:);
    else
        c = x;
    end
    %if any(c<0)
    %    error('all elements have to be greater than zero')
    %end
    
    % sort the vector in ascending order
    c = sort(abs(c),'ascend');
    c_norm = norm(c,1);
    N = length(c);
    sum = 0;
    
    for k = 1:N
        sum = sum + abs(c(k)) * ((N - k + 0.5)/N);
    end
    index = 1 - 2 * (sum/c_norm);
   
end
