% This code will plot the Time series, Average mutual information,
% Embedding dimension from FNN and Cao's method, 3-D phase portrait, 
% First return map, Recurrence plot and Degree distribution of a given signal
 
% We need to choose appropriate value of optimum time delay (from AMI) and
% minimum embedding dimension (either from FNN method or Cao's method)
 
% The source code files of ARFM_ami.m, ARFM_fnn, ARFM_cao.m, ARFM_delay_vec.m, ARFM_Recurrence_plot.m,
%  ARFM_Return_Map.m and ARFM_Visibility_Algorithm.m are required for the execution of the code
 
% Input:
% x:   Time sampled signal      [e.g., load('ABC.txt')]
% Fs:  Sampling Frequency (Hz)  (e.g., 10000)
% lag: Number of samples to be lagged while calculating AMI (Number in integer, e.g., 100)
% dmax: Maximum number of dimensions to be checked while calculating the embedding dimension (e.g., 20)
% Rtol: Distance threshold to decide false neighbor (Value lies between 10 < Rtol < 50, usually fixed to 10) (required for FNN code)
% Atol: Another criterion to remove false neighbors (Fixed value, 2)(required for FNN code) (Abarbanel et al. 1993) 
% threshold:    Threshold to detect recurrences (e.g. 10%, 20%, etc.)

% #### You can change the default values of the parameters according to the type of signal you studying ####

% References: 
% 1) Abarbanel, H.D., Brown, R., Sidorowich, J.J. and Tsimring, L.S., 1993. 
%    The analysis of observed chaotic data in physical systems. Reviews of modern physics, 65(4), p.1331.
 
% 2) Nayfeh, A.H. and Balachandran, B., 2008. Applied nonlinear dynamics: analytical, 
%    computational and experimental methods. John Wiley & Sons.

%% Setup
% Clear screen, close figures, clear variables from workspace
% clc; close all; clear all 
% Specify whether to truncate the data for speed or not
truncate = 'on'; % 'on' | 'off'

% Load the time series and truncate if N > 1000
% Load the time series
% x=load('ARFM_data/LCO.txt');
% x=load('ARFM_data/QP.txt');
% x=load('data_beta_4.mat');
% x=x.eta(1,:);

% x=u_f;
% t = t_GT_tau;

x=load('ARFM_data/Chaos.txt');
% x=load('ARFM_data/TAI.txt');
% x=load('ARFM_data/combustion_noise.txt');
% Input the sampling frequency
Fs=1/0.001;
% Calculate number of timesteps in the signal
N=length(x);
% Set the maximum number of timesteps in the signal
Nmax = 5000;
% Truncate the signal if it is too long (this can be omitted)
if (N>Nmax) && (strcmp(truncate,'on'))
    disp(['Reducing number of datapoints from ',num2str(N),' to ',num2str(Nmax),' in order to run quickly'])
    N=Nmax; x=x(end-N+1:end); t=t(end-N+1:end);
end

%% Plot the signal
% Centre the signal such that mean(x) = 0
x = x - mean(x);        
% Generate a time vector of the sampled data using a given value of sampling frequency (Fs)
delta_t = 1/Fs;   t = 1:N;  t = t*delta_t;

% Plot the signal
figure(201); %set(gcf,'Position',[450,450,950,500])
plot(t,x,'c','linewidth',1.5,'color',[.8 .1 0]);  set(gca,'fontsize',18)
title('The signal','Color', 'b','fontsize',20,'Fontweight','bold');
xlabel('Time (s)','fontsize',20)
ylabel('Amplitude','fontsize',20)
set(gca,'linewidth',1.5)
 
%% Plot average mutual information to get optimum time delay
% Set the maximum lag (in units of timesteps)
lag=100;
% Extracting the average mutual information and corresponding lag positions from the 'ARFM_ami.m' function
[v,L] = ARFM_ami(x,lag);
 % Lag is a non-dimensional time
% Convert lag values into milli-seconds
timelag = L/Fs*1e3;
 %% Plot the Average Mutual Information versus the time delay
[pks_AMI,tau_AMI]=findpeaks(-v);
Delay = tau_AMI(1);
tau = floor(Fs*Delay/1000);    % Division is by 1000 becasue the input is in 'ms'
disp(['The delay is ',num2str(tau_AMI(1)),' timesteps']) 
  %
figure(202); %set(gcf,'Position',[550,450,800,500])
plot(timelag, v, 'o-k','MarkerSize',8,'linewidth',1.5,'color',[0 .5 0]); %set(gca,'fontsize',18); grid on;
title('Average mutual information', 'Color', 'b')
xlabel('Lag (ms)')%,'fontsize',20)
ylabel('$\mathrm{I_{AB}}$')%,'fontsize',20)
set(gca,'linewidth',1.5)
%% Ask the user to input the optimum time delay
Delay = input('Input the optimum time delay (the first minimum of AMI) in ms:');
% convert to an integer number of timesteps
tau = floor(Fs*Delay/1000);    % Division is by 1000 becasue the input is in 'ms'
% Display to screen
% disp(['The delay is ',num2str(tau),' timesteps']) 

%% Obtaining embedding dimension from a False Nearest Neighbor method
% Set typical input values for algorithm
dmax=20; Rtol=10; Atol=2;
% Extract the percentage value of false nearest neighbors from ARFM_fnn.m
[FNN]=ARFM_fnn(x,t,Fs,dmax,tau,Rtol,Atol);
% plot the variation of percentage of false nearest neighbors with increase in embedding dimension
figure(203);
% set(gcf,'Position',[550,450,800,500])
plot(1:length(FNN),FNN,'-oy','MarkerSize',10,'linewidth',1.5,'color',[1 .5 0]); set(gca,'fontsize',20); grid on;
title('False nearest neighbors', 'Color', 'b','fontsize',20,'Fontweight','bold')
xlabel('Embedding dimension (d)','fontsize',22)
ylabel('False neighbors (%)','fontsize',22)
set(gca,'linewidth',1.5)

%% Obtaining embedding dimension from a Cao's method
[E1,E2] = ARFM_cao(x,tau,dmax); 
%% Plot for calculating embedding dimension from Cao's method
% In the plot quantities E1 & E2 are plotted against embedding dimension
figure(204); 
set(gcf,'Position',[550,450,800,500])
plot(E1,'-ok','MarkerSize',10,'linewidth',1.5); set(gca,'fontsize',18); grid on; hold on
plot(E2,'-rs','MarkerSize',10,'linewidth',1.5); hold off
title('Cao''s method', 'Color', 'b','fontsize',20,'Fontweight','bold')
xlabel('Dimension (d)','fontsize',20)
ylabel('E_{1} & E_{2}','fontsize',20)
set(gca,'linewidth',1.5)
%% Ask the user to provide an embedding dimension by visual inspection
d = input('Provide an appropriate embedding dimension from the comparison of FNN and Cao''s method and press Enter key ');
disp(['The embedding dimension is ',num2str(d)]) 

%% Phase space reconstruction
% Find length of the delay vectors (M) using the values of number of data points in the signal (N),
% minimum embedding dimension, (d) and the optimum time delay (tau) as M=N-(d-1)*tau
M=N-(d-1)*tau;
% Find number of delayed vectors of a signal from the given values of embedding dimension and delay
[Y]=ARFM_delay_vec(x,tau,d,M);  
 % 3-dimensional plot reconstructed phase portrait in embedding dimensional space
figure(205);
set(gcf,'Position',[550,450,800,500])
if size(Y,2) == 2
    plot(Y(:,1),Y(:,2),'k','linewidth',1.5,'color',[0 .5 0]);  set(gca,'fontsize',18); grid on;
    title('2D Phase Portrait', 'Color', 'b','fontsize',20,'Fontweight','bold')
else
    plot3(Y(:,1),Y(:,2),Y(:,3),'k','linewidth',1.5,'color',[0 .5 0]);  set(gca,'fontsize',18); grid on;
    title('3D Phase Portrait', 'Color', 'b','fontsize',20,'Fontweight','bold')
    zlabel('\itx\rm(t+2\tau)','fontsize',20)
end
xlabel('\itx\rm(t)','fontsize',20)
ylabel('\itx\rm(t+\tau)','fontsize',20)
set(gca,'linewidth',1.5)
drawnow

%% Plot First Return Map
pks= ARFM_return_map(x,Fs);

%% Plot Recurrence Plot
% Set the threshold (in %)
threshold=20;
S=ARFM_recurrence_plot(x,t,Fs,10,1,threshold);
 
%% Degree distribution of node from Visibility algorithm
deg=ARFM_visibility_algorithm(x,Fs);
