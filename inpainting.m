clc;
close all;
x = [-1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]; %Blog Data. -1 are liberals and 1 are conservatives
N = length(x);
x = transpose(x);
A = zeros(40); %Adjacency Matrix
A(1,3) = 1;% Adjacency matrix is custom made because for blog data the connections are made when one blog 
A(1, 10) = 1;% has links to other blogs. Therefore no automatic software can accomplish this task
A(3, 1) = 1;% However in case of audio signals the graph can be made using provided code
A(6, 10) = 1;
A(7, 40) = 1;
A(7, 15) = 1;
A(7, 10) = 1;
A(7, 20) = 1;
A(7, 19) = 1;
A(8, 25) = 1;
A(8, 9) = 1;
A(8, 17) = 1;
A(8, 19) = 1;
A(9, 16) = 1;
A(9, 8) = 1;
A(10, 8) = 1;
A(10, 9) = 1;
A(10, 6) = 1;
A(10, 19) = 1;
A(10, 7) = 1;
A(10, 15) = 1;
A(10, 17) = 1;
A(10, 20) = 1;
A(14, 25) = 1;
A(15, 16) = 1;
A(15, 10) = 1;
A(15, 8) = 1;
A(15, 7) = 1;
A(15, 20) = 1;
A(16, 14) = 1;
A(16, 15) = 1;
A(16, 9) = 1;
A(17, 8) = 1;
A(19, 17) = 1;
A(19, 10) = 1;
A(19, 8) = 1;
A(20, 10) = 1;
A(23, 32) = 1;
A(23, 28) = 1;
A(23, 35) = 1;
A(23, 31) = 1;
A(23, 26) = 1;
A(23, 36) = 1;
A(24, 31) = 1;
A(25, 14) = 1;
A(25, 26) = 1;
A(25, 8) = 1;
A(26, 21) = 1;
A(26, 25) = 1;
A(26, 23) = 1;
A(26, 24) = 1;
A(26, 38) = 1;
A(26, 31) = 1;
A(26, 40) = 1;
A(26, 29) = 1;
A(27, 38) = 1;
A(29, 36) = 1;
A(30, 32) = 1;
A(30, 31) = 1;
A(30, 33) = 1;
A(30, 37) = 1;
A(31, 26) = 1;
A(31, 23) = 1;
A(31, 28) = 1;
A(31, 30) = 1;
A(31, 33) = 1;
A(31, 32) = 1;
A(31, 38) = 1;
A(31, 37) = 1;
A(32, 26) = 1;
A(32, 31) = 1;
A(32, 33) = 1;
A(32, 30) = 1;
A(33, 26) = 1;
A(33, 39) = 1;
A(33, 38) = 1;
A(33, 31) = 1;
A(33, 30) = 1;
A(33, 36) = 1;
A(35, 39) = 1;
A(36, 28) = 1;
A(36, 33) = 1;
A(36, 23) = 1;
A(37, 31) = 1;
A(37, 30) = 1;
A(37, 38) = 1;
A(38, 33) = 1;
A(38, 26) = 1;
A(38, 36) = 1;
A(38, 39) = 1;
A(39, 35) = 1;
A(39, 33) = 1;
A(39, 38) = 1;
A(40, 7) = 1;
A(40, 26) = 1;

A = double(A);

outdegree = 0; %Assigning weight to each link. Equal weight is assigned using outdegree.
for i = 1:40
    for j = 1:40
        outdegree = outdegree + A(i, j); %Calculate outdegree. 
    end
    for k = 1:40
        if outdegree > 0
            A(i, k) = A(i, k)/outdegree; % Divide each link weight with outdegree to gice unifrom weight
        end
    end
    outdegree = 0;
end

xn =  [-1 0 -1 -1 0 -1 -1 -1 0 -1 -1 -1 -1 0 -1 -1 -1 -1 -1 0 1 0 1 1 0 1 1 1 0 1 1 1 1 0 1 1 1 1 1 0];
xn = transpose(xn); % This is the distorted signal

plot(x);
title('Original Signal');
grid on;
axis equal;
figure;
plot(xn);
title('Signal with unknown values');
grid on;
axis equal;

Z = ones(N, 2); % N = 40
X = zeros(N, 2); % These two vectors are being used to create a new vector which has all unknown values at
j = 1;% the end. And also these vectors will be used to rearrange the adjacency matrix
k = 1;

for i = 1:N
    if xn(i) == 0 % If the element in corrupted signal is zero
        Z(j, 1) = 0; % Then make the Z element 0 (the first column)
        Z(j, 2) = i; % And store the index of this elemeny in Z (second column)
        j = j + 1; %The index i of xn will help in rearranging the Adjacency matrix
    end
    if xn(i) ~= 0 % Else if the element is not zero
        X(k, 1) = xn(i); % Then do the same with the X vector
        X(k, 2) = i;
        k = k + 1;
    end
end

noz = j - 1; % Noz is the index upto which the Z vector will be sliced to get all the zero elements
ZA = Z(1:noz, 1); % Get all the zero elements of xn (using Z) and their indexes
XA = X(1:N-noz, 1); % Get all the one elements of xn (using X) and their indexes
index1 = Z(1:noz, 2); % The indexes were stored in the second columns of X and Z
index2 = X(1:N-noz, 2);

XN = vertcat(XA, ZA);% Vertically concatenate the XA and ZA to get sorted XN which has all the zero elements below the ones
index = vertcat(index2, index1); % Get their respective sorted indexes from xn

clear ZA XA index1 index2;

Y = zeros(N); % Using these two helper matrices new arranged matrix R will be created which will correspond
R = zeros(N); % to the signal XN which has all its zeros below the ones

for i = 1:N % N = 40
    j = index(i);
    Y(:, i) = A(:, j); % Copy columns according to the index vector which gives the index corresponding to 
end					% ones first and zeros second

for i = 1:N % Here the columns are in correct order but the rows are not acc. to the index vector.
    j = index(i);
    R(i, :) = Y(j, :); % Copy the rows into R matrix in arranged order
end

clear Y; % R is our new arranged matrix

I = eye(40);
Ai = ctranspose(I - R)*(I - R);

%Graph total variation regulation(10 unknown and 30 known)
Lamda = 0.01; 
L =  Lamda*Ai;
Im = eye(40); 
i=1;
j=1;
while i==j && i <= length(XN); %Making an identity matrix where the diagonal elements are zeros after the length of known part.
    if(i>30)
        Im(i,j) = 0;  
    end
    i = i+1;
    j = j+1;
end
phi = Im + L; 
xrecons = phi\XN; %Closed form solution of GTVR


for i = 1:length(xrecons) %Reconstructed signal values were threshold around zero. 
    if xrecons(i) >= 0 % So positive value is set to +1
        xrecons(i) = 1;
    end
    if xrecons(i) < 0 % Negative value is set to -1
        xrecons(i) = -1;
    end
end

H = zeros(1, length(xrecons));

for i = 1:length(xrecons)
    j = index(i);
    H(j) = xrecons(i);
end

%Plotting the recovered signal
figure;
plot(1:length(H) ,H);
title('GTVR')
grid on;
axis equal;
%Graph total variation minimization(10 unknown and 30 known);
Auu = Ai(N-noz+1:N, N-noz+1:N); %Separating unknown unknow matrix form Ai
AuuInv = inv(Auu); % Inverse of Unknown unknown matrix.
Aum = Ai(N-noz+1:N,1:N-noz); %Separating unknown know matrix from Ai
Xu = -(AuuInv*Aum*Xm); % Closed form of GTVM

for i = 1:length(Xu) %Reconstrued signal values are threshold around zeros
    if Xu(i) >= 0 % So positive value is set to +1
        Xu(i) = 1;
    end
    if Xu(i) < 0 %Negative value is set to -1
        Xu(i) = -1;
    end
end
figure;
plot(1:length(Xu) ,Xu);
title('GTVM (Only unknown values)')
grid on;
axis equal;