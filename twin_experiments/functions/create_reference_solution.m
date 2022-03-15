function [t_true, Aa_true] = create_reference_solution(P)
    rng('default')

    % Initial condition
    IC_eta  =   P.pert * ones(P.N_m,1);
    IC_etad =   IC_eta ;
    
    
%     IC_eta    =   zeros(P.N_m,1);
%     IC_eta(1) = 1;
%     IC_etad    =   zeros(P.N_m,1);
        
    IC_v    =   zeros(P.N_c,1);
    IC      =   [IC_eta; IC_etad; IC_v]; 
    t_v     =   (P.t_min:P.dt:P.t_max);
    % Add parameter if requested
    if isfield(P, 'param_est')
        if P.param_est == 1
            IC_param    =   P.beta;
            IC          =   vertcat(IC, IC_param);
        elseif P.param_est == 2
            IC_param    =   [P.beta; P.tau];
            IC          =   vertcat(IC, IC_param);
        end
    end

    [t_true, Aa_true] = ode45(@(t,y)gov_eqns_adv_PE(t,y,P), t_v, IC);
    % Force Y to be a matrix modes x time
    if size(Aa_true,1) > size(Aa_true,2)	
        Aa_true    =   Aa_true';    
    end
    
end