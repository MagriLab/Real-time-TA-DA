function [obs_idx, Filter] = fn_define_observation_indices(Filter, ESN_P)
    % ------------------------------------------------------------------- %
    t_start_i       =   Filter.t_start / Filter.dt_mic + 1;   
    t_stop_i        =   Filter.t_stop / Filter.dt_mic + 1;
    % ------------------------------------------------------------------- % 
    % Select when to start the washout
    t_start_ESN     =   Filter.t_start + ESN_P.N_washout * ESN_P.dt; 
    % Start time of state and parameter estimation 
    t_start_SPE      =   t_start_ESN;  
%     t_start_SPE      =   Filter.t_start + 0.1;  
    % Delay start time of parameter estimation only (i.e. drop SE)
    t_start_PE      =   t_start_SPE + 0.025; 
    
    % ------------------------------------------------------------------- %
    % Select assimilation frequencies
    kmeas_small     =   35; % (iniitlialisation)
    Filter.k_meas   =   35; %(general)
    % ------------------------------------------------------------------- %
    % Translate into indices
    t_start_SPE_i	=   round(t_start_SPE/ Filter.dt_mic) + 1;
    t_start_PE_i	=   round(t_start_PE/ Filter.dt_mic) + 1;
    t_stop_wash_i 	=   round(t_start_ESN/Filter.dt_mic) + 1;
    t_start_wash_i 	=   t_stop_wash_i - ...
                        double(ESN_P.N_washout) * ESN_P.dt/Filter.dt;
    % ------------------------------------------------------------------- %
%     tau_i       =   Filter.tau / Filter.dt;  
%     t_init      =   t_start_i-200*tau_i;
    t_init   	=   Filter.t_min  / Filter.dt_mic + 1;  
    part_0      =   [t_init,t_start_i];  
%     part_0      =   t_start_i-2*tau_i:1:t_start_i;    
    if t_stop_wash_i <= t_start_PE_i
        part_1	=   t_start_i:kmeas_small:t_start_wash_i;
        part_2	=   t_start_wash_i:1:t_stop_wash_i;
        part_3	=   t_stop_wash_i:kmeas_small:t_start_PE_i;
        part_4	=   t_start_PE_i:Filter.k_meas:t_stop_i;
    else 
        part_1	=   t_start_i:kmeas_small:t_start_PE_i;
        part_2 	=   t_start_PE_i:Filter.k_meas:t_start_wash_i;
        part_3 	=   t_start_wash_i:1:t_stop_wash_i;
        part_4 	=   t_stop_wash_i:Filter.k_meas:t_stop_i;
    end
    obs_idx     =   [part_0, part_1(2:end),part_2(2:end),part_3(2:end),part_4(2:end)];
    obs_idx     =   unique(obs_idx);
%     if length(obs_idx) ~= length(unique(obs_idx))
%         error('Check observation indices')
%     end
    % ------------------------------------------------------------------- %
    % Update Filter
    Filter.N_A              =   length(obs_idx);
    [~,Filter.init_SPE_i]	=   min(abs(t_start_SPE_i - obs_idx));
    [~,Filter.init_PE_i]    =   min(abs(t_start_PE_i - obs_idx));
    [~,Filter.init_BE_i]    =   min(abs(t_stop_wash_i - obs_idx));
end