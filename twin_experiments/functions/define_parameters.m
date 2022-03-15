function P = define_parameters(P)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Parameter definition
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isfield(P, 'assimilate_modes')
        if P.twin_experiment
            P.dt_exp = P.dt;
        end
        P.t1KF    	=   int64(P.t_minKF/P.dt_exp) + 1;   %if P.t1KF ==0; P.t1KF = 1; end
        P.tendKF   	=   int64(P.t_maxKF/P.dt_exp) + 1; 
        P.N_t   	=   (P.t_max - P.t_min) / P.dt ;    
        % Chack whether we are assimilating modes or microphone data
        if  P.assimilate_modes
            % Choose where to apply assimilation (1 - apply, 0 - not apply)
            P.y_0	=   [ones(1,P.N_m),...	% Observations for velocity modes
                        ones(1,P.N_m),...	% Observations for pressure modes
                        zeros(1,P.N_c)];  	% Observations for advection eq  
        else
            % Choose where to apply assimilation (1 - apply, 0 - not apply)
            P.y_0	=   [zeros(1,P.N_m),...	% Observations for velocity modes
                        zeros(1,P.N_m),...	% Observations for pressure modes
                        zeros(1,P.N_c)];  	% Observations for advection eq   
            if P.twin_experiment
                % Microphone parameters
                P.x_mic	=   linspace(0.25,.95,P.N_mic);
            end
        end
    else
        disp('No assimilate_modes field')
    end
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Initialise and store contant operand vectors and matrices 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    D_c         =   Cheb(P.N_c);          % Chebyshev matrix
    zeta        =   zeros(P.N_m,1) ;      % Damping factor
    cos_jpixf   =   zeros(1,P.N_m) ;      
    sin_jpixf   =   zeros(1,P.N_m) ;      
    jpi         =   zeros(P.N_m,1) ;      
    for j = 1:P.N_m
        zeta(j)         =   P.C1 * j^2 + P.C2 * sqrt(j) ;   
        jpi(j)          =   j * pi ;
        cos_jpixf(j)    =   cos(jpi(j) * P.x_f) ;  
        sin_jpixf(j)    =   sin(jpi(j) * P.x_f) ;  
    end
    % Store values  
    P.D_c       =   D_c;
    P.zeta      =   zeta ;      
    P.jpi       =   jpi;
    P.sin_jpixf =   sin_jpixf ;   
    P.cos_jpixf =   cos_jpixf ;
end