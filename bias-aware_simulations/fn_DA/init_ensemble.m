function [Af, py] = init_ensemble(p_mic, Filter)
% function [Aa, py, Filter] = init_ensemble(p_mic, Filter)
    % Initialise State
    [IC_ENS, py] = init_state(p_mic, Filter);
    % Initialise Parameters
    if Filter.E_Params
        % Add to the IC uniform random distribution of the beta and tau
        P_ENS   =   init_params(Filter);
        IC_ENS 	=   horzcat(IC_ENS, P_ENS);
    end
    Af     	=   IC_ENS;
end
% ======================================================================= %

function [IC_ENS, py] = init_state(p_mic, Filter)

    % Retreive variables:
    m       =   Filter.m;
    sig     =   Filter.sig;
    N_mic   =   Filter.N_mic;
    N_m     =   Filter.N_m;
    % Multivariate random distribution of initial pressure observation
    sig_f	=   p_mic * sig + 1E-4 ;
    Cpp_0   =   eye(N_mic).* mean(sig_f.^2);
    py      =   mvnrnd(p_mic, Cpp_0, m)'; 
    
    % Transform IC in pressure to IC in alphaj
    IC_alphaj   =   mvnrnd(ones(1,N_m), eye(N_m).* mean(sig_f.^2), m);
    
%     IC_alphaj   =   (- pinv(Filter.sin_omj_mic) * py)';
    
    IC_etaj     =   - IC_alphaj;     
%     IC_v        =   zeros(Filter.m, Filter.N_c);
    IC_v        =    zeros(Filter.m, Filter.N_c) + (Filter.cos_omjxf* IC_etaj')';

    IC_ENS      =   [IC_etaj, IC_alphaj, IC_v];
    
end
function p_ENS = init_params(P)
% Auxiliary function that initialises the parameters to estimate as an
% uniform random distribution. 
% ----------------------------------------------------------------------- %
    IC_k    =   P.sig_PE;
    if P.E_Params == 3            
        p_ENS = [unifrnd(P.beta * (1-IC_k), P.beta * (1+IC_k), P.m, 1),...
                unifrnd(P.tau * (1-IC_k), P.tau * (1+IC_k), P.m, 1),...
                unifrnd(P.C1 * (1-IC_k), P.C1 *(1+IC_k),P.m,1),...
                unifrnd(P.C2 * (1-IC_k), P.C2 * (1+IC_k), P.m, 1)];
    elseif P.E_Params == 2            
        p_ENS = [unifrnd(P.beta * (1-IC_k), P.beta * (1+IC_k), P.m, 1),...
                unifrnd(P.tau * (1-IC_k), P.tau * (1+IC_k), P.m, 1)];
    elseif P.E_Params == 1
        p_ENS = unifrnd(P.beta * (1-IC_k), P.beta * (1+IC_k), P.m, 1);
    end
% ----------------------------------------------------------------------- %
end