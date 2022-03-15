function [pks]= ARFM_return_map(x,Fs)
% This code will help in exploring the properties of the higher dimensional
% system into a lower dimensional subspace.
% This method is also known as Poincaré section 
% It helps in revealing the various types of motion such as periodic, 
% quasiperiodic and chaotic exhibited by the given dynamical system 
% In the case of periodic systems of known period T, the first return map 
% would produce N discrete points
% in the two-dimensional phase space. Here, N could be an integer with a 
% value equal to 1 for periodic
% oscillations, 2 for period-2 oscillations, and so on.
% If N is a non-integer but rational number, the motion is called 
% mode-locked, whereas, if the N is irrational, then the motion is called 
% quasiperiodic. For the mode-locked state, the first return map will show 
% a set of discrete points in the intersection plane. For quasiperiodic 
% motion, the first return map shows a closed dense circle in the 
% two-dimensional subspace. The Poincaré sections of chaotic or 
% higher-dimensional quasiperiodic systems do not exhibit a simple
% geometric pattern, as observed for the periodic or 2-period quasiperiodic 
% motions . In order to create first return map, 
%       (1) the extreme events (such as local maxima or minima or zero
%            crossing) of every cycle of the signal is calculated and 
%       (2) then the first event is plotted against the next one. 
% ______________ Input _______________________
% x:    Time series
% Fs:   Sampling frequency 
% ____________Reference ______________________
% Nayfeh, A.H., and Balachandran, B., 2008. Applied nonlinear dynamics: analytical,
%       computational and experimental methods. John Wiley & Sons.
 
% E.g.: [pks]= ARFM_Return_Map(x,10000)
%% Get first return map
 
% Calculate all positive maxima's of every cycle of the given time series
% Here, pks - values of local maxima's and locs - time indices corresponding to these maximum values
[pks,locs]=findpeaks(x,'MinPeakDistance',10,'Minpeakheight',0);
 
figure(305);
set(gcf,'Position',[550,450,800,500])
% set(gcf,'Position',[300,300,1400,500]);
% subplot (121)
 
% Create a reference diagonal line in the plot with size sufficiently higher than maximum value of the signal
M=max(pks)+2*max(pks);
plot(0:1/Fs:M,0:1/Fs:M,'linewidth',1.5,'Color', [0 .5 0]);  grid on;  hold on;
 
% Plot a first return map wherein the first maxima of the signal is plotted with the next maxima
plot(pks(1:end-1),pks(2:end),'.k','Color', [.8 .1 0],'MarkerSize',14); set(gca,'fontsize',18)
 
axis tight
title('First return map','fontsize',22,'Color', 'b','Fontweight','bold')
xlabel('\itx\rm_{max} (i)','fontsize',24)
ylabel('\itx\rm_{max} (i+1)','fontsize',24)
xlim([0,M]);
ylim([0,M]);
set(gca,'linewidth',1.5)
