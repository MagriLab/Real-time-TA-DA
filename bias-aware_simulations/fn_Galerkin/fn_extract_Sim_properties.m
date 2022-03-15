function P = fn_extract_Sim_properties(Sim, P)
    % =================================================================== %
    u1_0    =   Sim.Mean.u1;
    u2_0    =   Sim.Mean.u2;    
    p2_0    =   Sim.Mean.p2;    
    T1_0	=   Sim.Mean.T1;
    T2_0    =   Sim.Mean.T2;
    % Geometry properties
    P.x_f  	=   Sim.Geom.Lu;                % Flame location        [m]
    P.L_0  	=   Sim.Geom.Lb + Sim.Geom.Lu;  % Rijke tube length     [m]
    % Compute the weighted averaged temperature and mean flow velocity
    x       =   P.x_f / P.L_0;
    T_0     =   (1 - x) * T2_0 + x * T1_0;	% MF temperature        [K]
    u_0     =   (1 - x) * u2_0 + x * u1_0;  % MF velocity           [m/s]
    p_0     =   p2_0;                       % MF pressure           [Pa]
    % Compute the other thermodynamic properties from p and T
    gamma	=   Sim.Mean.gamma;            	% Heat capacity ratio   [-]
    R       =   8.314 / 0.02895;            % Gas constant
    c_0     =   sqrt(gamma * R * T_0);      % MF speed of sound     [m/s]
    rho_0   =   p_0 / (R * T_0);            % MF density            [kg/m3]
    % Store mean flow values
    P.Mean.T_0      =   T_0;
    P.Mean.u_0      =   u_0;
    P.Mean.p_0      =   p_0;
    P.Mean.c_0      =   c_0;
    P.Mean.rho_0	=   rho_0;
    P.Mean.gamma	=   gamma;
    % =================================================================== %
    % Heat release rate properties
    P.law   =   Sim.Mean.Heat_law;
    P.beta  =   Sim.Mean.beta;
    P.tau   =   Sim.Mean.Tau;
    P.kappa =   Sim.Mean.kappa;
    % =================================================================== %
    P.t_min     =   0;                      % Minimum integration time  [s]                          
    P.t_max     =   Sim.Measurement.Time;	% Maximun integration time  [s]                    
    P.Fs_mic    =   Sim.Measurement.Fs;     % Measurement timestep  [s]
    P.dt        =   Sim.dt;                 % Integration timestep      [s]
    % =================================================================== %
    P.x_mic     =   Sim.Measurement.Mic_Pos + P.x_f; % Mics location    [m]
    P.N_mic     =   length(P.x_mic);
end

