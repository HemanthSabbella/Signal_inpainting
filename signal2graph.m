function adjacency = signal2graph(im1)
% Takes discret signal as input
% Gives Adjacency matrix as output

[m,n] = size(im1); 
                                         % m is row of the given signal
                                         % n is the clounms of the given sinal 

adjacency = ones(m*n,m*n);
                                         % Creates a matrix with ones of size m*n x m*n

for i = 1:m*n
    x1=im1(1,i);
    for j = i+1:m*n
        x2=im1(1,j);
        adjacency(i,j) = exp((x1-x2)^2); % Makes the upper tringle of adjacency matrix
        adjacency(j,i) = exp((x2-x1)^2); % Makes the lower tringle of adjacency matrix
    end
end

                                         % We dont need to calculate diagonal values as they are always 
