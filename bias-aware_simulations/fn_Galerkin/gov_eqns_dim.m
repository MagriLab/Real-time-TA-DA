function dYdt = gov_eqns_dim(t, y, P, law )

    % Get parameters
    N_m         =   P.N_m ;
    N_c         =   P.N_c ;
    
    zeta        =   P.zeta ;
    D_c         =   P.D_c ;
    u_0         =   P.Mean.u_0;
    rho_0     	=   P.Mean.rho_0;
    c_0         =   P.Mean.c_0;
    
    gamma      	=   P.Mean.gamma;
    L_0         =   P.L_0;
    omega_j     =   P.omega_j;
    cos_omjxf   =   P.cos_omjxf ;
    sin_omjxf   =   P.sin_omjxf ;
    
    % Rename variables
    ETAj  	=   y(1 : N_m);
    ALPHAj	=   y(N_m+1 : 2*N_m);
    if ~isfield(P,'E_Params') || ~P.E_Params % P.param_est == 0 
        K_p     =   [];
        beta   	=   P.beta ;
        tau    	=   P.tau ;
    elseif P.E_Params == 1
        beta    =   y(end);
        tau    	=   P.tau ;
        K_p     =   0;
    elseif P.E_Params == 2
        beta    =   y(end - 1);
        tau     =   y(end);
        K_p     =   [0;0];
    end
    
    
    % ================= Apply BCs for the advection eqns ================ %   
    V    	=   y(2*N_m+1 : 2*N_m+N_c);
    u_f     =   cos_omjxf * ETAj;       % Velocity at the flame at t
    V2      =   [u_f; V];
    u_ftau  =   V2(end);                % Velocity at the flame at t - tau
    if u_ftau > 1E10
        wtf=1;
    end
    if length(y) >= 2*N_m+2*N_c
        W    	=   y(2*N_m+N_c+1 : 2*N_m+2*N_c);
        p_f     =   - sin_omjxf * ALPHAj;   % Pressure at the flame at t
        W2      =   [p_f; W];
        p_ftau  =   W2(end);                % Pressure at the flame at t - tau
    end 
    

    % ===================== Define heat release law ===================== %    
    if contains(law,'Q13') ||strcmp(law,'sqrt') || isempty(law)
        qdot    =	beta * (sqrt(abs(u_0/3 + u_ftau)) - sqrt(u_0/3)) ; %[kg/s3] = [W/m2]
        qdot    =   - 2 * qdot * (gamma - 1) / L_0 * sin_omjxf'; % [Pa/s]
        
    elseif contains(law,'sqrt_mod')
        qdot    =	beta * (sqrt(abs(u_0/3 + u_ftau*3)) - sqrt(u_0/3)) ; %[kg/s3] = [W/m2]
        qdot    =   - 2 * qdot * (gamma - 1) / L_0 * sin_omjxf'; % [Pa/s]
        
    elseif contains(law,'Q1')
        qdot	=   beta * (sqrt(abs(1 + u_ftau)) - sqrt(1)) ;
        qdot    =   - 2 * qdot * (gamma - 1) / L_0 * sin_omjxf';
    elseif contains(law,'tan')
        kappa   =   P.kappa;
        qdot    =   beta * sqrt(beta/kappa) * atan(sqrt(kappa/beta) .* u_ftau); % [m/s3]
%         qdot   =   - (beta^2 * p_ftau) / (beta + kappa * u_ftau^2) * ...
%                 cos_omjxf / sin_omjxf; % [Pa/s]
            
        qdot    =   - 2 * qdot * (gamma - 1) / L_0 * sin_omjxf';
    elseif contains(law, 'cub')
        kappa   =   P.kappa;
        qdot   =   - (beta^2 * p_ftau) / (beta + kappa * u_ftau^2) * ...
                cos_omjxf / sin_omjxf; % [Pa/s]
    else
        error(['Heat release law ''',law,'''not defined'])
    end
    
    % ================= Dimensional acoustic equations ================== %
    K_ETAj   	=   (omega_j ./ (rho_0 * c_0)).* ALPHAj;
    K_ALPHAj	=   - (omega_j * rho_0 * c_0) .* ETAj ...
                    - c_0/L_0 * zeta.* ALPHAj ...
                    + qdot ;   
                
    K_V        	=     2 / tau * D_c * V2 ; % Chebyshev
    if length(y) >= 2*N_m+2*N_c
        K_W        	=     2 / tau * D_c * W2 ; % 
    end
    if any(isnan(K_V))
        wtf=1;
    end

    % ========================== Output vector ========================== %
    if length(y) >= 2*N_m+2*N_c
        dYdt = [K_ETAj; K_ALPHAj; K_V(2:end); K_W(2:end); K_p]; 
    else
        dYdt = [K_ETAj; K_ALPHAj; K_V(2:end); K_p]; 
    end
end
