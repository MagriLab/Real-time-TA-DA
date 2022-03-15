function [E1,E2] = ARFM_cao(x,tau,d)
% This code will help in finding minimum embedding dimension using a method proposed by Cao (1997)
% _____________ Inputs ______________________
% x:    The time series
% tau:  Optimum time delay obtained from average mutual information (AMI)
% d:    Maximum number of embedding dimensions required for the calculation
% _____________ Reference__________________ 
% Reference:
% Cao, L., 1997. Practical method for determining the minimum embedding dimension of a scalar time series.
% Physica D: Nonlinear Phenomena, 110(1-2), pp.43-50.
% ________________________________________ 
% We need to get a plot of variation of parameters E1(d) and E2(d) (defined in the paper by Cao, 1997) with dimension 
% This plot will give the appropriate value of minimum embedding dimension required for phase space reconstruction 
% In this case, E1(d) will stop changing once all the false neighbors are resolved after a dimension d_0. 
% Then, d_0+1 can be chosen as the minimum embedding dimension required for the construction of the phase space.
%________________________________________
% Initializing the null vectors of E1 and E2 for d+1 dimension
E_m = zeros(d+1,1);
Estar_m = zeros(d+1,1);
% Calculating the values of E1 and E2 for every embedding dimension
N = length(x);
Nmax = 10000;
if (N>Nmax) %&& (strcmp(truncate,'on'))
    disp(['Reducing number of datapoints from ',...
        num2str(N),' to ',num2str(Nmax),' in order to run quickly'])
    N=Nmax; x=x(end-N+1:end);% t=t(end-N+1:end);
end

for m = 1:d+1
    % Reduce the data vector to the nearest multiple of delay
    n = floor(length(x)/tau)*tau;    
    % Store the delayed signals i.e., Y(d)={x(t), x(t+tau)... } in the form of row vectors,
    % e.g., for d=2 you have x(t) and x(t+tau), and data length is [n-2*tau]    
    Y = zeros(n-m*tau, m);  % Initializing a delay vector, where n-m*tau is the total number of delay vectors
    % Loop to calculate delay vectors for a given time series using Takens delay embedding theorem (1981)
    for i = 1:m+1       
        Y(:,i) = x((i-1)*tau+1:(n-(m+1-i)*tau));
    end
    
    N = size(Y,1); % Length of the delay vectors
    % Initializing two vectors required to store the values of a(i,d) (Eq. (1) in Cao, 1997) 
    % and distance vector of Eq. (4) in Cao, 1997
    dist = zeros(N,1);
    a2 = zeros(N,1);
    % Loop to calculate the value of a(i, d) for every d as given by equation (1) in Cao (1997)
    for i = 1:N
        
        temp = Y(i,1:end);  % Swapping of Y vector after each iteration
        Y(i,1:end) = Y(1,1:end);
        Y(1,1:end) = temp;    
        % Obtain a distance matrix of points on the phase space trajectory
        % i.e., calculating the distance of every initial reference point on the phase space trajectory 
        % with respect to other points of trajectory and repeating the same process for N-1 points
        Rd = transpose(sqrt(sum((repmat(Y(1,1:end-1),size(Y,1)-1,1)...
            - Y(2:end,1:end-1)).^2,2)));   
        % Finding the nearest neighbor to a given reference point of the phase space trajectory i.e., finding the 
        % value of minimum distance from the distance matrix and the position of this neighbor from the reference point        
        [val_mRd, pos_mRd] = min(Rd);          
        % Obtain the distance nearest neighbor with other points on the trajectory due to increasing the dimension by 1.
        % This is needed for the calculation of E2        
        dist(i) = abs(Y(1,end) - Y(1+pos_mRd,end));
        % Avoid division by 0
        val_mRd = val_mRd + eps;        
        % Obtain the total distance on increasing the dimension by 1
        minRdplus1 = sqrt(val_mRd*val_mRd + dist(i)*dist(i)) + eps;               
        % Take the ratio of distances. This is needed for E1
        a2(i) = minRdplus1/val_mRd;     % Equation (1) in the paper by Cao (1997)
    end
       E_m(m) = mean(a2);          % Equation (2) in the paper
       Estar_m(m) = mean(dist);    % Equation (4) in the paper by Cao (1997)
end
      % Finding E1(d) and E2(d) for every dimension
      E1 = E_m(2:end)./E_m(1:end-1);              % Equation (3) in the paper by Cao (1997)
      E2 = Estar_m(2:end)./Estar_m(1:end-1);      % Equation (5) in the paper by Cao (1997)
end
