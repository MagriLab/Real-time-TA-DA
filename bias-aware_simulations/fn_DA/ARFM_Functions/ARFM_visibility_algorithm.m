function [deg,fig]=ARFM_visibility_algorithm(x,t,Fs)
% This code extracts the complex network for any given signal using a methode of "Visibility Graph Algorithm" (Lacasa et al. 2008)
% This code will help in plotting the degree distribution and getting an adjacency matrix required for 
% generating a complex network (using Gephi software) of any given time series 
% _____________ Inputs _____________________ 
% x:    Time series
% Fs:   Sampling frequency 
% e.g. [deg]=ARFM_Visibility_Algorithm(x,Fs); 
% ___________  Reference __________________
% Lacasa, L., et al. "From time series to complex networks: The visibility graph." 
%          Proceedings of the National Academy of Sciences 105.13 (2008): 4972-4975.
%%  
N=length(x);            % Calculate length of the signal 

% Obtain vector of local maxima of a given signal. Each maximum is 
% considered as a 'node' in the complex network and the line connecting 
% each node is called 'connection'. Find magnitude and location of all the 
%local peaks of a given signal    
[pks,~] = findpeaks(x,'MinPeakHeight',0,'MinPeakDistance',20);   
peak=pks';         % Obtain a vector of nodes from the signal 
% Obtain adjacency matrix. It is required to store the information about the connectivity between the nodes (Donner et al. 2010)
% If two nodes are connected then it is represented by '1' in the adjacency matrix and otherwise, it is referred as '0'
% Two peaks are connected if they are visible    ### For 'visibility condition' read: L. Lacasa et al. 2008,
% The two nodes are connected only if "there exists a straight line
% that connects these nodes and this line is not intersected by any of the intermediate points"
 
L1=length(peak);    % Find total number of nodes in the signal
 
A=zeros(L1);        % Defining a adjacency matrix
 
for i=1:L1-1
    for j=i+1:L1
        % Condition of visibility..... [Eq. (1) in Lacasa et al (2008)]
        if peak(i+1:j-1)<(peak(j)+(peak(i)-peak(j))*(j-(i+1:j-1))/(j-i))    
            A(i,j)=1;   % If the visibility condition holds, represent this connection by 1; otherwise, it will be automatically taken as 0
                        % To avoid self-connection A(i,i) is always considered as zero
        end
    end
    A(i,i+1)=1;         % Condition representing the connection between neighboring nodes
end
 
% The final adjacency matrix will be an upper triangular matrix
A=A+A';     % Generate a symmetric matrix of A
 
%%  Analyzing complex networks using Gephi software
 
% In order to analyze signal using Gephi software, first we need to create
% a matrix whose first row and the first column will be numbered from
% 1,2,3,..... to the length of the adjacency matrix. 
%This can be achieved by shifting the adjacency matrix in both row as well as column by one 
 
B(2:length(A)+1,2:length(A)+1)=A;         % Shift the Adjacency matrix by one row and one column
B(1,2:length(A)+1)=1:length(A);           % Number the first column, other than first element, as 1 to length of the Adjacency matrix
B(2:length(A)+1,1)=1:length(A);           % Number the first row, other than first element, as 1 to length of the Adjacency matrix
csvwrite('D.csv',B);                      % This creates an excel file containing information of connection of all the nodes of a signal. 
                                          % Use this file directly in the Gephi software to create a complex network
 
%% Find the degree distribution of a complex network mapped from the time series
 
% Degree of a node: sum of all the nodes that are connected to a given node
deg=sum(A);   % Sum all the numbers of adjacency matrix in column wise i.e., the sum all the number of connections to a particular node
 
table=tabulate(deg(2:L1-1)); % creates a frequency table of the data present in 'deg' vector
% This command arranges the data in a table as:
% 1st column ‚Äî a vector containing unique values (no repetition) of a 'deg' vector arranged in ascending order
% 2nd column ‚Äî The number of instances values in 1st column are repeated in a 'deg' vector
% 3rd column ‚Äî The percentage of each value in the total set of degree of a node
% http://www.mathworks.com/access/helpdesk/help/toolbox/stats/tabulate.html

% assigning each column of a table by different variable names
value=table(1:end-1,1);
count=table(:,2);
Percent=table(1:end-1,3);
 
x=find(Percent);            % Find the indices of a percentage vector having nonzero elements
value=value(x);             % Get elements from a 'value' vector having non-zero percentage
Percent=Percent(x);         % Get elements from a 'percentage' vector having non-zero percentage
 
% Plot degree distribution of complex network
% This is log-log (to the base 10) plot between the percentage of nodes with k number of
% connections and the degree of the node
 
logvalue = log10(value);
logPercent = log10(Percent);
 
fig = figure(7);
plot(logvalue,logPercent,'ok','MarkerSize',8)
 
LV=logvalue(2:end); LP=logPercent(2:end);
C=polyfit(LV,LP,1);  % Linear fit to the degree distribution plot
Lin_fit=C(1)*LV+C(2);
hold on;
% C(1)
plot(LV,Lin_fit,'r-.','MarkerSize',10,'linewidth',1.5)
title('Degree distribution');
xlabel('$\log(k)$','FontSize',20)
ylabel('$\log(P(k))$','FontSize',20)
text(max(LV-.2),max(Lin_fit)-.1,sprintf('Slope = %0.2f',C(1)), 'horizontalAlignment', 'center')
% set(gca,'Position',[0.6 .25 0.3 0.6])
