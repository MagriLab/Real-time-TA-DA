function [t, Af] = fn_forecast_ESN(t, Af, dt_f, Filter)


    % ============= Define function of governing equations ============== %
    fun	=   @(t,psi) gov_eqns_dim(t,psi,Filter,Filter.law);
    

    % =========== Propagate each ensemble member individually =========== %
    if size(Af,1) == Filter.m
        for j = 1:size(Af,1)
            psi_j       =   Af(j,:)';
            [~, psi_j] 	=   RK4_solver(t,psi_j,dt_f,fun);
            Af(j,:)   	=   psi_j;
        end
        t   =   t + dt_f;
    else
        [t, Af] 	=   RK4_solver(t,Af,dt_f,fun);
    end
end