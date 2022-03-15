function dYdt = gov_eqns_adv_PE(t, y, P )

    % Get parameters
    N_m         =   P.N_m ;
    N_c         =   P.N_c ;
    cos_jpixf   =   P.cos_jpixf ;
    sin_jpixf   =   P.sin_jpixf ;
    jpi         =   P.jpi ;
    zeta        =   P.zeta ;
    D_c         =   P.D_c ;
    beta        =   P.beta ;
    tau         =   P.tau ;
    
    % Rename variables
    ETAj        =   y(1 : N_m);
    ETADj_jpi   =   y(N_m + 1 : 2 * N_m);
    V           =   y(2 * N_m + 1 : 2 * N_m + N_c);
    
    if ~isfield(P,'param_est') || P.param_est == 0 
        K_p     =   [];
    elseif P.param_est == 1
        beta    =   y(end);
        K_p     =   0;
    elseif P.param_est == 2
        beta    =   y(end - 1);
        tau     =   y(end);
        K_p     =   [0;0];
    end
    
%     disp(['beta = ', num2str(beta), 'tau = ', num2str(tau)])
    
    u_f     =   cos_jpixf * ETAj;   % Velocity at the flame at t
    V2      =   [u_f; V];
    u_ftau  =   V2(end);            % Velocity at the flame at t - tau
    
    % Square heat release law
    Q       =     beta * (sqrt(abs(1/3 + u_ftau)) - sqrt(1/3)) ;
%     Q       =     beta * (sqrt(abs(1 + u_ftau)) - sqrt(1)) ;
    
    % Galerkin discretisation 
    K_ETAj          =   jpi.* ETADj_jpi;
    K_ETADj_jpi     =   - jpi.* ETAj - zeta.* ETADj_jpi - 2 * Q * sin_jpixf';
    
    % Chebyshev discretisation
    K_V             =     2 / tau * D_c * V2 ;

    dYdt = [K_ETAj; K_ETADj_jpi; K_V(2:end); K_p]; 
end
