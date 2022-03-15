function Sim = fn_analysis_wave(i, Sim)
    % This version includes analysis update to the history
    Filter = Sim.Filter;    
    % ====================== INITIALISE ENSEMBLE ======================== %
    if i == 1
        fprintf('Initialising ensemble'); 
        Sim = init_ensemble(Sim);
        fprintf('\nAnalysis progress: 0%% '); 
    % ==================== ASSIMILATE OBSERVATIONS ====================== %
    elseif i > 1
        if mod(i,round(Filter.num_analysis/10))==0
            fprintf([num2str(ceil(i/Filter.num_analysis*100)),'%% '])
        end
    else
        if mod(i,round(Filter.num_analysis/10))==0
            fprintf([num2str(ceil(i/Filter.num_analysis*100)),'%% '])
        end
        tA  =   Filter.analysis_index(i);
        pA  =   zeros(Filter.m,Filter.N_mic); 
        for mi = 1:Filter.m
            [~,pA(mi,:)]  =   fn_p_from_GH(mi, Sim);
        end
        % ---------------- Create state matrix with past ---------------- %
        Mean    =   Sim.Mean;
        J       =   max([Mean.Tau, Mean.Tau+Mean.Tu, Mean.Tau+Mean.Td]);
        J       =   ceil(J / Sim.dt);
        
        % Select 10 points of the history
        N_pts   =   11;
        t_hist  =   linspace(Sim.ts(tA-J),Sim.ts(tA),N_pts);
        Gs_hist =   interp1(Sim.ts,Sim.Ensemble.Gs', t_hist)';
        Hs_hist =   interp1(Sim.ts,Sim.Ensemble.Hs', t_hist)';
        
        Af      =   [Gs_hist, Hs_hist];
        
%         Af  =   [Sim.Ensemble.Gs(:,tA),...
%                  Sim.Ensemble.Hs(:,tA)];
             
        % ---------- Augment matrix with parameters and pressure -------- %
        if Filter.param_estimation
            Af  =   [Af,...
                     Sim.Ensemble.Taus(:,tA)];
        end
        Af  =   [Af,pA]';
        % -------------------- Retrieve observations -------------------- %
        tA_exp  =   Filter.analysis_index_exp(i);
        py      =   Sim.Truth.P_mic(:,tA_exp);
        % ------------------------ Appply EnSRKF ------------------------ %
        [Aa_m, Psi_acomb, py] = EnSRKF(Af, py, Filter);
        Aa_new  =   Aa_m + Psi_acomb;
        
        % ------------- Store analysis state and pressures -------------- %
        Sim.Ensemble.Gs_phy(:,tA)   =   Aa_new(1,:)';
        Sim.Ensemble.Hs_phy(:,tA)   =   Aa_new(N_pts+1,:)';  
        
        ts_hist     =   Sim.ts(tA-J:tA);
        Gs_hist_new	=   spline(t_hist,Aa_new(1:N_pts,:)',ts_hist);
        Hs_hist_new	=   spline(t_hist,Aa_new(N_pts+1:2*N_pts,:)',ts_hist);
        
        Sim.Ensemble.Gs(:,tA-J:tA)	=   Gs_hist_new;
        Sim.Ensemble.Hs(:,tA-J:tA) 	=   Hs_hist_new;    
        
        if Filter.param_estimation
            Sim.Ensemble.Taus(:,tA)	=   Aa_new(2*N_pts+1,:);
        end
        Filter.Observations{i} = py;
        Filter.pA_f{i}      =   pA; 
        Filter.pA_a{i}      =   Aa_new(end-Filter.N_mic+1:end,:)';
    end
    
    Sim.Filter =    Filter;
end

% =========================================================================
function [Aa_m, Psi_acomb, py] = EnSRKF(Af_, py, Filter)
% Algorithm for the Ensemble Square-Root Kalman filter. The inputs required
% are the forecast state matrix (Af_), the observations (py), the
% observations mapping matrix (M), the observations covariance matrix (Cee)
% and the struct of parameter of the problem (P).
% -------------------------------------------------------------------------
    m       =   Filter.m;
    N_mic	=   Filter.N_mic;
    % Create observations with some deviation from the truth
    psi_f_m     =   mean(Af_, 2); 
    sig_mic     =   py * Filter.sig_mic + 1E-5;
    Cee         =   eye(Filter.N_mic).* mean(sig_mic.^2); 
    py          =   mvnrnd(py, Cee);
    % Observable quantities
    y_0     =   zeros(size(Af_,1),1);
    y_0(end-N_mic+1:end) = 1; 
    q = N_mic;
    % Observation mapping matrix
    M     	=   zeros(q,length(y_0));
    iq    	=   1;
    for i = 1:length(y_0)
        if y_0(i) ~= 0  
            M(iq,i) =   1;      
            iq      =   iq + 1;    
        end
    end    
    % Create matrices of observations and mean
    Y  	=   py';
    Af_m	=   psi_f_m;         
    % Compute the analysis state mean value
    Psi_f   =   (Af_ - Af_m);
    if ~isreal(Psi_f)
        Psi_f   =   Af_ - Af_m + 1E-5;
        disp('Not real')
    end
    % Transform into observation space
    S       =   M * Psi_f; 
    W       =   (m - 1) * Cee + S * S.'; 
    % Check condition number
    try cond_W  =   cond(W);
        if cond_W > 1e15
            K       =   S.' * pinv(W);
        else
            K       =   S.'/W;
        end
    catch
        disp('ILL CONDITIONED!')
    end
    % Compute the analysis state mean value
    Aa_m    =   Af_m +  Psi_f * K * (Y - M * Af_m);
    % Obtain analysis ensemble perturbations
    VEV         =   K * S;
    [V,E]       =   eig(VEV,'matrix');   V  =   real(V); E   =   real(E) ;
    Psi_acomb	=	Psi_f * V * ((eye(m)  - E)^(0.5)) * V.';
end

% =========================================================================
function Sim = init_ensemble(Sim)
% TODO - Initialise as function of pressure and not Gs Hs                   % TO DO
    Filter  =   Sim.Filter;
    Truth   =   Sim.Truth; 
    Mean    =   Truth.Mean;
    Geom    =   Truth.Geom;
    
    % Retreive variables:
    m       =   Filter.m;
    sig     =   Filter.sig;
    Tau     =   Mean.Tau;
    Tu      =   Mean.Tu;
    dt      =   Truth.ts(2);
    N       =   Geom.N;
    
    % Time indices for the truth
    tf_f   	=   Filter.analysis_index_f(1);
    tf_s 	=   Filter.analysis_index_s(1);
    t0_f    =   1;
    t0_s    =   tf_s - tf_f + 1;    
    
    len     =   length(t0_f:tf_f);
    
    % True values
    G      =   Truth.Gs(t0_s:tf_s)';
    H      =   Truth.Hs(t0_s:tf_s)';
    X      =   Truth.Xis(:,t0_s:tf_s);
    Q     =   Truth.Qs(t0_s:tf_s)'; 

%     % mvnrnd OPTION  ---------------------------------------- %
%     vec 	=	vertcat(G, H, Xis, Qs);
%     sig     = 	vec * Filter.sig + 1E-4 ;
%     var     =   mean(sig.^2,2);
%     
%     IC_ENS  =   mvnrnd(vec, eye(size(vec)) .* sig, m);
%     all_vec =   repmat(reshape(vec,[1,size(vec)]), m,1,1) ;
%     dIC     =   all_vec - IC_ENS;
%     IC_ENS  =   all_vec + dIC;
%                                 
%     Sim.Ensemble.Gs(:,t0:tf) = squeeze(IC_ENS(:,1,:));
%     Sim.Ensemble.Hs(:,t0:tf) = squeeze(IC_ENS(:,2,:));
%     Sim.Ensemble.Qs(:,t0:tf) = repmat(Qs,[m,1]);
%     Sim.Ensemble.Xis(:,:,t0:tf) = repmat(Xis_re,[m,1,1]);
%     % -------------------------------------------------------- %
    
    
    % UNIFRN OPTION --------------------------------------------%
%     IC_G    =   sort([G(end) * (1-sig),G(end) * (1+sig)]);
%     IC_G    =   unifrnd(IC_G(1), IC_G(2), m,1);
%     IC_H    =   sort([H(end) * (1-sig),H(end) * (1+sig)]);
%     IC_H    =   unifrnd(IC_H(1), IC_H(2), m,1);
%     IC_Q    =   sort([Q(end) * (1-sig),Q(end) * (1+sig)]);
%     IC_Q    =   unifrnd(IC_Q(1), IC_Q(2), m,1);
%     Sim.Ensemble.Gs(:,t0_f:tf_f) = G - (IC_G - G(end));
%     Sim.Ensemble.Hs(:,t0_f:tf_f) = H - (IC_H - H(end));
%     Sim.Ensemble.Qs(:,t0_f:tf_f) = Q - (IC_Q - Q(end));
    
    % CORRELATED UNIFRN
    IC_fact     =   unifrnd(1-sig, 1+sig, m,1);
    
    
    for j=1:N
        Sim.Ensemble.Xis(:,j,t0_f:tf_f) = X(j,:) - X(j,end) * (IC_fact - 1);
    end
    
    Sim.Ensemble.Gs(:,t0_f:tf_f) = G - G(end) * (IC_fact - 1);
    Sim.Ensemble.Hs(:,t0_f:tf_f) = H - H(end) * (IC_fact - 1);
    Sim.Ensemble.Qs(:,t0_f:tf_f) = Q - Q(end) * (IC_fact - 1);
    
%     % fOR NOW IM KEEPING QS AND XIS EQUAL TO THE TRUTH AS I AM NOT UPDATING  %% TO DO
%     % THEM WITH THE FILTER
%     Xis      =   reshape(Truth.Xis(:,t0_s:tf_s), [1,N,len]);
%     Sim.Ensemble.Qs(:,t0_f:tf_f)    =     repmat(Qs,[m,1]);
%     Sim.Ensemble.Xis(:,:,t0_f:tf_f) =     repmat(Xis,[m,1,1]);
    % -------------------------------------------------------- %
    
    if isfield(Sim.Ensemble,'Taus')
        Tau     =   Mean.Tau;
        IC_k    =   Filter.sig_PE;
        IC_Tau  =   unifrnd(Tau * (1-IC_k), Tau * (1+IC_k), Filter.m, 1);
        Sim.Ensemble.Taus(:,t0_f:tf_f) = repmat(IC_Tau, [1,len]);
    end
    
end