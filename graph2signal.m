function signal = graph2signal(adjacency,value,point)
% Takes ajacency matrix, one known value, and its index in the signal as input
% Gives discret signal as output

[m,n] = size(adjacency);
% m is row of the given adjacency matrix
% n is the clounms of the given adjacency matrix

signal = zeros(1,n);
                                 % Creates a matrix with zeros of size m*n x m*nxn

for i = 1:n
    x = log(adjacency(point,i)); % takes log of that particular value from matrix
    y = sqrt(x);                 % Takes the square root of the above value
    signal(1,i) = value - y;     % Substract the above value from the gven known value to obtain the orignal signal
end


