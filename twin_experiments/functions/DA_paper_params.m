function param  =   DA_paper_params(fig)
    % This function defines the parameters required to replicate the
    % figures in Nóvoa, A. & Magri, L. (2021). Real-time thermoacoustic 
    % data assimilation.
    param = [];
    param = base_params(param);    
    if fig ==  5 || fig ==  6  % Modes no PE timeseries or QP error
        param.assimilate_modes	=   true;   % Assimilate the modes?
        param.param_est     	=   0;      % Estimate beta and/or tau? 
    elseif fig == 7 % Modes PE params
        param.assimilate_modes	=   true;   
        param.param_est     	=   2;      
    elseif fig == 9 || fig ==  10 % Mics no PE timeseries or QP error
        param.assimilate_modes	=   false;  
        param.param_est     	=   0;      
    elseif fig == 11
        param.assimilate_modes	=   false;  
        param.param_est     	=   2;      
        param.k_meas            =   1.5;
    elseif fig == 12
        param.assimilate_modes	=   false;  % Assimilate the modes?
        param.param_est     	=   2;
        param.k_meas            =   1;
        param.N_mic             =   15;
    elseif any(fig == [13,14,15]) % Chaos butterfly timeseries and error
        param.assimilate_modes	=   true;
        param.param_est     	=   0;
    elseif fig == 16 || fig == 17 % Chaos timeseries or PSD
        param.param_est     	=   0;      
        param.m             	=   100;  	
        param.k_meas          	=   0.5; 
    elseif fig == 18 % Chaos PE params with 8 microphones
        param.param_est     	=   2;       
        param.m             	=   300;  	
        param.k_meas          	=   .5; 
        param.inflation         =   1.2;       
    elseif fig == 19 % Chaos PE params with 15 microphones
        param.param_est     	=   2;      
        param.m             	=   300;  
        param.k_meas          	=   .5; 
        param.inflation         =   1.02;    
        param.N_mic             =   15;
    end
    
    param.k_meas 	=   param.k_meas / param.dt ;

end

function P = base_params(P)
    % Thermoacoustic model parameters
    P.C1           	=   0.1 ;   % Damping parameter C1 = 0.1 
    P.C2           	=   0.06 ;  % Damping parameter C2 = 0.06   
    P.tau          	=   0.2;    % Time delay
    P.x_f          	=   0.2;    % Flame location
    
    % Numerical model parameters
    P.N_m          	=   10;     % Number of Galerkin modes
    P.N_c          	=   10;     % Number of Chebishev modes    
    P.N             =   P.N_m * 2 + P.N_c;     	% State vector size
    P.dt           	=   0.001;  % Time step
    P.pert         	=   0.005;  % Initial condition for Galerkin modes
    P.t_min        	=   0;      % Initial integration time    
    P.m            	=   10;     % Ensemble size
    
    % EnSRKF parameters
    P.twin_experiment   =   true;   % Twin experiment?
    P.IC_param_fact     =   0.25;   % Initial uncertainty for parameters
    P.sig               =   0.25;   % Initial uncertainty on the state 
    P.k_meas            =   2;      % Time between analysis
    P.N_mic             =   6;      % Number of microphones
    P.mic_sig_frac      =   0.01;   % Microphone uncertainty
    P.inflation         =   1.0;    % Inflation parameter
end