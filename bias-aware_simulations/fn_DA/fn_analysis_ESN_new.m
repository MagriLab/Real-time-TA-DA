function [Aa, py, Filter] = fn_analysis_ESN_new(i, Af, pA, p_mic, Filter)

    % ====================== INITIALISE ENSEMBLE ======================== %
    if i == 1
        [Af_augmented, py] = init_ensemble(p_mic, Filter);
        fprintf('\nAnalysis progress: '); 
    end
    % ==================== ASSIMILATE OBSERVATIONS ====================== %
    % ---------------------- Retrieve observations ---------------------- %
    sig_f	=   p_mic * Filter.sig_mic/10 + 1E-4 ;
    Cee     =   eye(Filter.N_mic).* mean(sig_f.^2);
    py      =   mvnrnd(p_mic, Cee);      
    % -------------------------- Appply EnSRKF -------------------------- %
    [Aa_m, Psi_acomb, py] = EnSRKF(Af_augmented, py, Cee, Filter);
    Aa  =   Aa_m + Psi_acomb;  
    % ---------------- Inflate if unphyscal parameters ------------------ % 
    if Filter.E_Params
        if i > 15
            condition = check_bounds(Aa, Filter);        
            if any(condition)
                rho         =   Filter.inflation;
                Af_m        =   mean(Af_augmented, 2);
                Af_inflated	=   Af_m + rho * (Af_augmented - Af_m);
                Aa        	=   Af_inflated;
            end
        else
            Param_i         =   Filter.N-Filter.E_Params+1:Filter.N;
            Aa(Param_i,:)	=   Af_augmented(Param_i,:);
        end
    end
    % ----------------------------- Output ------------------------------ % 
    if any(i == 1:round(Filter.N_A/5):Filter.N_A)	% Print progress
        fprintf([' ',num2str(round(i/Filter.N_A*100)),'%%'])
    end
     % Force output to be (m x N) and remove pressure
    if size(Aa, 1) ~= Filter.m
        Aa = Aa';
    end
    Aa  =   Aa(:,1:end-Filter.N_mic);
end


%% ===================================================================== %%

function condition = check_bounds(Aa,P)
% Auxiliary function that defines the allowed values for the estimated
% parameters. 
% ----------------------------------------------------------------------- %
    condition   =   false;
    if P.E_Params == 3
        condition   =   [min(Aa(P.N-3,:))<-1;...
                        max(Aa(P.N-3,:)) > 10; ...
                        min(Aa(P.N-2,:))<-1;...
                        max(Aa(P.N-2,:)) > 1; ...
                        min(Aa(P.N-1,:))<0.2;...
                        max(Aa(P.N-1,:)) > 10; ...
                        min(Aa(P.N,:)) < 0.005; ...
                        max(Aa(P.N,:)) > 0.8,...
                        ];
    elseif P.E_Params == 2
        condition   =   [min(Aa(P.N-1,:))<10;... % beta
                        max(Aa(P.N-1,:)) > 1E15; ...
                        min(Aa(P.N,:)) < 0.0001; ...
                        max(Aa(P.N,:)) > 0.2];
    elseif P.E_Params == 1
        condition   =   [min(Aa(P.N,:))<0.05; ...
                        max(Aa(P.N,:)) > 10];
    end
end