
function [t_comp, Aa_comp] = create_comparison_solution(P, Aa)

    %%%%%%%%%%%%%%%%%%% Solution without Data Assimliation %%%%%%%%%%%%%%%%%%%%
    % Initial condition
    IC_comp     =	Aa(1,:,1);
    t_v         =   (P.t_minKF:P.dt:P.t_max);

    [t_comp, Aa_comp] = ode45(@(t,y)gov_eqns_adv_PE(t,y,P), t_v, IC_comp);
    % Force Y to be a matrix modes x time
    if size(Aa_comp,1) > size(Aa_comp,2)	
        Aa_comp    =   Aa_comp';    
    end
end