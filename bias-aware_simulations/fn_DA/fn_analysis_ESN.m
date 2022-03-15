function [Aa, py, Filter] = fn_analysis_ESN(i, Af_augment, p_mic, Filter)

    m           =   Filter.m;
    N           =   Filter.N;
    N_m         =   Filter.N_m;
    N_c         =   Filter.N_c;
    N_mic       =   Filter.N_mic;
    Ni_params  	=   2*N_m+N_c+1:N;
    % ====================== INITIALISE ENSEMBLE ======================== %
    if i == 1
        [Af, py]	=   init_ensemble(p_mic, Filter);
        pA       	=   - Filter.sin_omj_mic * Af(:,N_m+1:2*N_m)';
        % Only assimilate state
        [Af_augment, Filter] = fn_create_augmented_Af(Af',pA,1,0,Filter);
        Aa  =   Af_augment;
        py  =   mean(py,2);
    else
    % ==================== ASSIMILATE OBSERVATIONS ====================== %
    % ---------------------- Retrieve observations ---------------------- %
        sig_f	=   p_mic * Filter.sig_mic/10 + 1E-4 ;
        Cee     =   eye(N_mic).* 1E-10;%mean(sig_f.^2);
        py      =   mvnrnd(p_mic, Cee);      
        % -------------------------- Appply EnSRKF -------------------------- %
        [Aa_m, Psi_acomb, py] = EnSRKF(Af_augment, py, Cee, Filter);
        Aa  =   Aa_m + Psi_acomb;  
        % ---------------- Inflate if unphyscal parameters ------------------ % 
        if Filter.E_Params
            if i > Filter.init_SPE_i
                condition       =   check_bounds(Aa, Filter);        
                if any(condition)
                    rho         =   Filter.inflation;
                    Af_m        =   mean(Af_augment, 2);
                    noise       =   mvnrnd(Af_augment'*0, diag(abs(Af_m*1E-5)))';
                    Af_inflated	=   Af_m + rho * (Af_augment - Af_m)  + noise;
                    N_last_p	=   size(Aa,1) - Filter.N_mic;
                    N_first_p	=   N_last_p - Filter.E_Params + 1;

                    spread      =   std(Af_inflated(N_first_p:N_last_p,:),1,2)...
                                  ./min(Af_inflated(N_first_p:N_last_p,:),[],2);

                    condition  	=   check_bounds(Af_inflated, Filter);     
                    if any(condition)
                        wtf=1;
                    end
                    if any(spread > 0.4) || any(condition)
                        Aa      =   Af_augment;
                    else
                        Aa    	=   Af_inflated;
                    end
                end
            end
        end
    end
    % ----------------------------- Output ------------------------------ % 
     % Force output to be (m x N) and remove pressure
    if size(Aa, 1) ~= m
        Aa = Aa';
    end
    Aa  =   Aa(:,1:end-N_mic);
    % If first iteration and parameter estimation, return also parameters
    % initial values
    if i == 1 && Filter.E_Params
        Aa      =   horzcat(Aa, Af(:,Ni_params));
    end
end


%% ===================================================================== %%

function condition = check_bounds(Aa,P)
% Auxiliary function that defines the allowed values for the estimated
% parameters. 
% ----------------------------------------------------------------------- %
    condition   =   false;
    N_last_param      =   size(Aa,1) - P.N_mic;
    if P.E_Params == 3
        condition   =   [min(Aa(P.N-3,:))<-1;...
                        max(Aa(P.N-3,:)) > 10; ...
                        min(Aa(P.N-2,:))<-1;...
                        max(Aa(P.N-2,:)) > 1; ...
                        min(Aa(P.N-1,:))<0.2;...
                        max(Aa(P.N-1,:)) > 10; ...
                        min(Aa(P.N,:)) < 0.005; ...
                        max(Aa(P.N,:)) > 0.8...
                        ];
    elseif P.E_Params == 2
        condition   =   [min(Aa(N_last_param-1,:))<P.beta/10;... % beta
                        max(Aa(N_last_param-1,:)) > P.beta*10; ...
                        min(Aa(N_last_param,:)) < P.tau/5; ... % tau
                        max(Aa(N_last_param,:)) > P.tau*1.2];
    elseif P.E_Params == 1
        condition   =   [min(Aa(N_last_param,:))<P.beta/5; ... % beta
                        max(Aa(N_last_param,:)) > P.beta*5];
    end
end