function [t_KF, Aa_KF, Cpp_tr] =   EnSRKF_parfor(P, Aa_true, t_true)

    rng('default')  % For reproducibility

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1. Create the true reference solution. Compute pressure at the
    % microphone locations if using microphones
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    if ~P.assimilate_modes
        sin_mic     =   zeros(P.N_mic,P.N_m);
        for mici = 1:P.N_mic
            sin_mic(mici,:)    =   sin(P.jpi.' .* P.x_mic(mici));
        end
        p_mic_true	=  - sin_mic * Aa_true(P.N_m+1:2*P.N_m,:);
        pinv_mics   =   pinv(sin_mic);
    end
    k_meas  =   P.k_meas;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 2. Initialise the KF parameters and matrices
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    k_t     =   (P.t_max - P.t_minKF)/P.dt + 1;	% Length time space
    % Create the ensembles. Pressure states at different locations created 
    % from the referrence state using a multivariate distribution 
    if P.assimilate_modes
        sig         =   mean(abs(Aa_true(1:2*P.N_m,P.t1KF))) * P.sig; % Standard dev
        Cpp_0       =   eye(P.N_m) * sig^2;                         % Cov matrix for IC
        eta_ENS     =   mvnrnd(Aa_true(1:P.N_m,P.t1KF),Cpp_0,P.m);  
        etad_ENS 	=   mvnrnd(Aa_true(P.N_m+1:2*P.N_m,P.t1KF),Cpp_0,P.m);  
    else
        sig_forecast	=   p_mic_true(:,P.t1KF) * P.sig + 1E-4 ;
        Cpp_0           =   eye(P.N_mic).* mean(sig_forecast.^2);
        IC_p            =   mvnrnd(p_mic_true(:,P.t1KF),Cpp_0,P.m)'; 
        etad_ENS        =   (- pinv_mics * IC_p)';
        % Use the same magnitude but different sign for velocity modes
        eta_ENS     =   - etad_ENS;
            % Initialise advection modes to zero
        %     IC_pmtau    =   mvnrnd(p_mic_true(:,int64(P.t_minKF-P.tau/P.dt + 1)),Cpp_0,P.m)'; 
        %     eta_dot_true	=   (- pinv_mics * IC_pmtau)';
        %     uf              =   P.cos_jpixf * eta_ENS;
        %     uftau           =   P.cos_jpixf * - eta_dot_true;
        %     v_ENS         	=   init_advection(uf, uftau, P.N_c, P.m);
    end
    v_ENS           =   zeros(P.m, P.N_c);
    IC_ENS      	=   [eta_ENS, etad_ENS, v_ENS];
    % Modify parameters when applying parameter estimation
    if P.param_est   
        IC_k    =   P.IC_param_fact;
        P.N     =   P.N + P.param_est;   % State vector w/ parameter estimation
        % Add to the IC uniform random distribution of the beta and tau
        p_ENS	=   [unifrnd(P.beta * (1-IC_k), P.beta *(1+IC_k),P.m,1),...
                   	unifrnd(P.tau * (1-IC_k), P.tau * (1+IC_k), P.m, 1)];
        IC_ENS	=   horzcat(IC_ENS, p_ENS(:,1:P.param_est));
        P.y_0	=   horzcat(P.y_0, zeros(1,P.param_est));
    end
    if ~P.assimilate_modes
        % Modify y_0 to account for microphones
        P.y_0       =   horzcat(P.y_0, ones(1,P.N_mic));
    end
    q   	=   sum(P.y_0~=0);      % Number of observations    
    % Initialise matrices for EnSRKF
%     Af   	=   zeros(P.N,P.m);             % Forecast ensmeble modes matrix
    Af_m	=   zeros(length(P.y_0),P.m);  	% Pressure forcast mean matrix
    M     	=   zeros(q,length(P.y_0));   % Observation operation matrix 
    iq    	=   1;
    for i = 1:length(P.y_0)
        if P.y_0(i) ~= 0  
            M(iq,i) =   1;      
            iq      =   iq + 1;    
        end
    end
    Y       =   zeros(q,P.m);     	% Matrix of identical observ. columns
    try
        Aa      =   zeros(k_t,P.N,P.m);	% Analysis Solution matrix
    catch
        pack
        Aa      =   zeros(k_t,P.N,P.m);	% Analysis Solution matrix
    end
    TX      =   zeros(k_t,1);      	% Time vector
    Cpp_tr  =   zeros(1,k_t);       % Trace of the Ensemble covariance
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 3. EnSRKF loop
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Aa_t        =   IC_ENS(:,:)';   % Initial consition for integration
    Cpp_tr(1)   =   trace(Cpp_0);
    ti_next     =   1;
    time_vec_KF         =   P.t_minKF:P.dt:P.t_maxKF;
    num_analysis_steps  =	floor(length(time_vec_KF) / P.k_meas);
    for ii = 1:num_analysis_steps-1
        % Update timesteps
        ti	=   ti_next;    % Current initial time, base 1
        if  ii == 1
            ti_next     =   P.k_meas; % Current final time, base 1
            ti_exp      =   P.t1KF + k_meas;   % Exp time for assimilation, base t1KF
        else
            ti_next     =   ti + P.k_meas;
            ti_exp      =   ti_exp + k_meas;
        end
        t_vec_ode   =   time_vec_KF(ti:ti_next);  % Integration time
        % Predict state at next measurement timestep
        Az_par = [];
        parfor mi = 1:P.m
            [~, Az]	=	ode45(@(t,yy)gov_eqns_adv_PE(t,yy,P),...
                            t_vec_ode, Aa_t(:,mi));
            Az_par(:,:,mi)  =   Az;
        end
        
        Af                  =   squeeze(Az_par(end,:,:)); 
        Aa(ti:ti_next,:,:)	=   Az_par ;    
        TX(ti:ti_next)    	=   t_vec_ode;
        
        % Ensemble covariance, mean and forecast ensemble perturbations
        Af_modes        =   Af(1:2*P.N_m,:);
        Cpp             =   cov(Af_modes.');
        Cpp_tr(ti_next)	=   trace(Cpp);
 

        if P.assimilate_modes
            % Compute the pressure at the different locations and add to state vec
            Af_     	=   Af;
            psi_f_m     =   mean(Af_, 2); 
            sig         =   mean(abs(Aa_true(1:2*P.N_m,ti_exp)),1) * P.sig;
            Cee         =   eye(q) .* sig^2;
            py          =   mvnrnd(Aa_true(1:2*P.N_m,ti_exp),Cee); 
        else
            % Compute the pressure at the different locations and add to state vec
            Af_     	=   vertcat(Af, - sin_mic * Af(P.N_m+1:2*P.N_m,:));
            psi_f_m     =   mean(Af_, 2); 
            % Update the microphone covariance
            sig_mic     =   (p_mic_true(:,ti_exp)) * P.mic_sig_frac + 1E-5;
            Cee         =   eye(P.N_mic).* mean(sig_mic.^2); %diag(mean(sig_mic.^2));
            % Obtain the pressure observations gaussian from reference state
            py          =   mvnrnd(p_mic_true(:,ti_exp), Cee);
        end
        % Create matrices of observations and mean
        for yi = 1:P.m
            Y(:,yi)  	=   py;
            Af_m(:,yi)	=   psi_f_m; 
        end
        % Compute the analysis state mean value
        Psi_f   =   (Af_ - Af_m);
%         Psi_f   =   P.inflation * (Af_ - Af_m);
        if ~isreal(Psi_f)
            Psi_f   =   Af_ - Af_m + 1E-5;
            disp('Not real')
        end
        % Transform into observation space
        S       =   M * Psi_f; 
        W       =   (P.m - 1) * Cee + S * S.'; % ADD NOISE - SMALL MATRIX AS COMPARED TO THE TWO TERMS
        % Check condition number
        try cond_W  =   cond(W);
            if cond_W > 1e15
                K       =   S.' * pinv(W);
            else
                K       =   S.'/W;
            end
        catch
            disp('ILL CONDITIONED!')
            disp(P.m)
            disp(Aa_t)
        end

        % Compute the analysis state mean value
        Aa_m    =   Af_m +  Psi_f * K * (Y - M * Af_m);
        % Obtain analysis ensemble perturbations
        VEV         =   K * S;
        [V,E]       =   eig(VEV,'matrix');   V  =   real(V);    E   =   real(E) ;
        Psi_acomb	=	Psi_f * V * ((eye(P.m)  - E)^(0.5)) * V.';
        % Compute the analysis state (mic pressure) and transform into modes
        Aa_new      =   Aa_m + Psi_acomb;
        % Correct and store for next IC
        condition   =   false;
        % SKIP ANLYSIS IF PARAMS OUT OF BOUNDS-----------------------------
        if P.param_est == 2
            if P.beta == 0.2
                lb = 0.05;
            else
                lb = 0.2;
            end
            condition   =   [min(Aa_new(P.N-1,:))<lb;...
                            max(Aa_new(P.N-1,:)) > 2*P.beta; ...
                            min(Aa_new(P.N,:)) < 0.005; ...
                            max(Aa_new(P.N,:)) > 0.5];
        elseif P.param_est == 1
            condition   =   [min(Aa_new(P.N,:))<lb; ...
                            max(Aa_new(P.N,:)) > 10];
        end
        if any(condition)
            Aa(ti_next,:,:)	=   Af_m(1:P.N,:) + P.inflation * Psi_f(1:P.N,:);
            Aa_t(:,:)       =   Aa(ti_next,:,:);
    %         create_plot_per_member(Aa_true,param,time_vec_KF, ti_next)
        else
            Aa(ti_next,:,:)	=   Aa_new(1:P.N,:);
            Aa_t(:,:)       =   Aa(ti_next,:,:);
        end
    end %while loop % -----------------------------------------------------


    
    if ti_next < length(time_vec_KF)
        t_vec_extra     =   time_vec_KF(ti_next):P.dt:t_true(end);
        IC_extra        =   mean(Aa_t,2);     
        % Predict state at next measurement timestep
        [t_noDA, Aa_noDA] = ode45(@(t,y)gov_eqns_adv_PE(t,y,P),...
                            t_vec_extra, IC_extra);
                        
        Aa(ti_next:end,:,:)    =   repmat(Aa_noDA, [1,1,P.m]) ;    
        TX(ti_next:end)        =   t_noDA;
    end

    %% Store output  vectors
        Aa_KF   =   Aa;
        t_KF    =   TX;
end