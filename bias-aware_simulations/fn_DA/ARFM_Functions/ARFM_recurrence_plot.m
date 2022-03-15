function [S,t1, t2, eps]=ARFM_recurrence_plot(X,T,d,tau,threshold)
% This code will display recurrence plot of a given time series 
% In order to get recurrence plot, we first need to reconstruct a phase space using appropriate
% values of time delay and embedding dimension. After the reconstruction of phase space, the distance
% every point of the trajectory is calculated with respect to other points to create a distance matrix.
% Using an appropriate value of threshold (percentage of the maximum size of the attractor), a recurrence
% matrix is created. A two-dimensional plot of the recurrence matrix is called recurrence plot. 
% It is better to use less data point in order to reduce computation time (e.g. N=1000 points) 
% ________________ Inputs ____________________
% x:    Time series
% Fs:   Sampling frequency
% d:    Minimum embedding dimension
% tau:  Optimum time delay (in samples)
% threshold:    Threshold to detect recurrences (e.g. 10%, 20%, etc.) 
% _______________ Output ________________
% S:    Distance matrix
% e.g. [S]=Recurrence_plot(x,10000,6,20,20); 
% _______________ Reference _____________
% N. Marwan, M. Carmenromano, M. Thiel, and J. Kurths,
%       "Recurrence plots for the analysis of complex systems", Phys. Rep. 438, 237 (2007).
 
%% Plot the signal
N=20 / (T(2) - T(1));            % Calculate length of the signal
x = X(1:10:N);       
t = T(1:10:N);

N=ceil(length(x));
% dt = t(2)-t(1);
% tau = (tau*dt)+t(1);

%% Construction of delay vectors using Taken's embedding theorem
 
% Find length of the delay vectors (M) using known values of number of data 
% points in the signal (N), minimum embedding dimension (d), and optimum 
% time delay (tau) as: M=N-(d-1)*tau
 
M=N-(d-1)*tau;
 
% Initializing a matrix of delay vectors
Y = zeros(M,d);
 
% Loop to get time delayed vectors using a time delay embedding theorem of Taken's
for i=1:d
    Y(:,d-i+1) = x(1+(i-1)*tau:M+(i-1)*tau);
end
 

%% Get distance matrix
 
L1 = length(Y(:,1));   % Length of the column delay vecor
 
% Creating a distance matrix by calculating the distance of every point of phase space trajectory
% with other points, and repeating the same process for all other points
 
x1 = repmat(Y, L1, 1); % Replicating the M1 matrix L1 times
 
x2 = reshape(repmat(Y(:), 1, L1)', L1 * L1, d);  % Arranging delay vector matrix in such a way that every row of the matrix is repeated for L1 times
 
 
S = sqrt(sum( (x1 - x2) .^ 2, 2 ));         % Subtraction of x1 and x2 will create a distance vector
S = reshape(S, L1, L1);                     % Reshaping the distance vector into a distance matrix form

%% Get recurrence plot
 
% Creating the time vector of length L1 (i.e., number of points in the column of a distance matrix)
t1 = t(1:L1);
t2 = t1;
 
% Recurrence threshold can be given in terms of percentage of maximum size of the attractor (i.e., maximum value of distance from the distance matrix)
% We can use another threshold such as the percentage of Euclidian distance or fixed number of nearest neighbors etc. [Ref: Marwan (2007)]
eps=threshold/100;   % Set given recurrence threshold in percentage
 


