function [t_KF, Aa_KF, Obs, t_bias, bias, to, U_open] =   EnSRKF_Galerkin(Truth, Filter, ESN_P)

    % Print type of Assimilation

    % Retrieve true variables
    t_true  =   Truth.t_mic;
    p_true  =   Truth.p_mic;
    % Retrieve filter variables
    t_max   =   Filter.t_max;
    t_stop  =   Filter.t_stop;
    m       =   Filter.m;
    N       =   Filter.N;
    N_c     =   Filter.N_c;
    N_m     =   Filter.N_m;
    N_mic   =   Filter.N_mic;
    dt      =   Filter.dt;
    dt_f    =   dt;
    N_A   	=   Filter.N_A;
    
    sin_omj =   Filter.sin_omj_mic;
    
    % Retrieve ESN vsriables
    dt_esn      =   ESN_P.dt;
    N_washout   =   ESN_P.N_washout;   

    %
    % Initialise output matrices 
    [ti, ti_mic]	=   deal(1);
    [t, t_mic]   	=   deal(t_true(1));
    p_mic         	=   p_true(ti_mic,:); 
    
    t_KF    =   t;
    Aa_KF   =   zeros(m, N, 1);
    Af      =   zeros(m, N);
    Obs     =   zeros(N_mic, 1);
    t_Obs   =   t_mic;
    t_bias  =   [];
    
    U       =   zeros(1, N_mic);
    bias    =   zeros(N_mic, 1);
    
    % Define index location of alpha, psi and parameters
    Ni_alpha  	=   N_m+1:2*N_m;
    Ni_state	=   1:2*N_m+N_c;
    Ni_params  	=   2*N_m+N_c+1:N;
    indices     =   {Ni_state, Ni_params};
    % Select a number of stes to perform state estimation only
    SE        	=   0; 
    PE        	=   0;
    fprintf('\nAnalysis progress: '); 
    while true
        if any(ti_mic == 1:round(N_A/5):N_A)	% Print progress
            fprintf([' ',num2str(round(ti_mic/N_A*100)),'%%'])
        end
        % ==================== Loop over observations =================== % 
        % -------------- Compute forecast pressure states --------------- %
        pA    	=   - sin_omj * Af(:,N_m+1:2*N_m)';
%         pA_m  	=   mean(pA,2);
%         pA   	=   (pA - pA_m) + (pA_m + U); % Deviations + updated mean
        pA   	=   pA + bias(:,end); % Deviations + updated mean
        if sum(sum(isnan(pA)))
            wtf=1; % Some pressure calculation gives NaN!
        end
        [Af_aug, Filter]	=   fn_create_augmented_Af(Af,pA,SE,PE,Filter);
        % ------------ Perform batch initialisation/analysis ------------ %
        if any([SE,PE]) || ti_mic == 1 
            [Aa, py, Filter] = fn_analysis_ESN(ti_mic, Af_aug, p_mic, Filter);
            if isreal(Aa)
                % ----- Update state and/or parameters with analysis ---- %
                if ti_mic == 1 
                    Aa_KF(:,:,ti)	=   Aa;
                elseif ti_mic <= Filter.init_SPE_i % SE only for initialising state
                    Aa_KF(:,Ni_state,ti)	=   Aa(:,Ni_state);
                elseif ti_mic <= Filter.init_PE_i % SE&PE if PE for initialising params
                    Aa_KF(:,:,ti)	=   Aa;           
                else % Remove SE if not required
                    keep_is	=   cell2mat(indices(find([SE,PE])));
                    Aa_KF(:,keep_is,ti)	=   Aa;
                end
            else
                wtf=1; % Analysis matrix is not real
            end 
        else
            Cee     =   eye(N_mic).* 1E-10;
            py      =   mvnrnd(p_mic, Cee);     
        end
        if ti_mic == Filter.init_SPE_i
            PE	=   Filter.E_Params;
            SE	=   1;
        elseif  ti_mic == Filter.init_PE_i
            SE	=   Filter.E_State;
        end 
        Af	=   Aa_KF(:,:,ti); 
        % --------------------- Store observations ---------------------- %
        Obs(:,ti_mic)	=   py;
        t_Obs(ti_mic)   =   t_mic;
        % get next time for analysis and break if no more data is available
        try 
            ti_mic  =   ti_mic + 1;  
            t_mic   =   t_true(ti_mic);
            p_mic   =   p_true(ti_mic,:); 
        catch
            break
        end         
        % =================== Forecast bias with ESN ==================== % 
        if ti_mic >= Filter.init_BE_i + 1
            if ti_mic == Filter.init_BE_i + 1
                    % ----------- Initialise ESN in open loop ----------- %
                    p_hf	=   - sin_omj * squeeze(mean(Aa_KF(:,Ni_alpha,:),1));
                    [U_wash, to]   =   create_washout(dt_esn, 2*N_washout, ...
                                            t_KF, p_hf, t_Obs, Obs);
                    [U_open, ra]	=   advance_ESN('open', ESN_P, U_wash);
                    U               =   U_open(end,:)';
%                     % %%%%%%%%%%%% PLOT DEBUG
%                         figure; hold on
%                         l=plot(t_w,U_wash(:,1),'k-','linewidth',3);l.Color(4)=.2;
%                         plot(t_w, U_open(:,1), 'k*--')
%                     % %%%%%%%%%%%% PLOT DEBUG
                % --------------- Store first bias values --------------- %
                bias(:,1)	=   U;
                t_bias(1) 	=   t_Obs(end); 
                t2_esn      =   1;
            else
            % ------------ Update bias with observatons data ------------ %
                pA_m	=   - sin_omj * mean(Aa_KF(:,Ni_alpha,end),1)';
                U       =   py' - pA_m;
                bias(:, end) 	=   U;
            end
            % --------------- Forecast ESN in closed loop --------------- %
            t_esn  	=   t_bias(end):dt_esn:t_mic;% From current time to next DA
            Nti     =   length(t_esn);
            [U, ra] =   advance_ESN('closed', ESN_P, ra, U, Nti);
%                 %%%%%%%%%%%% PLOT DEBUG
%                     plot(t_esn, squeeze(U(1,:)), '*-')
%                     plot(t_Obs(end), Obs(1,end) - p_hf(1,end), 'kx','MarkerSize',12)
%                 %%%%%%%%%%%% PLOT DEBUG
            % ------------------------ Store bias ----------------------- %
            t1_esn  =   t2_esn + 1;
            t2_esn  =   t1_esn + Nti - 2;
            bias(:, t1_esn:t2_esn) 	=   U(:,2:end);
            t_bias(t1_esn:t2_esn)   =   t_esn(2:end);            
        end
        % ================ Loop over forecast time steps ================ % 
        while true
            % ------------ Forecast members individually ------------ % 
            [t, Af]     =   fn_forecast_ESN(t, Af, dt, Filter);
            ti          =   ti + 1;     
            % Store
            Aa_KF(:,:,ti)	=   Af;
            t_KF(ti)        =   t;
            if	abs(t - t_mic) < 1E-10
                break
            end
        end
    end
    % ============== Integrate further as mean if required ============== %    
    if t_stop < t_max
        Af  =   mean(Af,1)';
        tj  =   0;
        Aa_extend	=   zeros(1, N, 1);
        t_extend    =   zeros();
        while t < t_max
            [t, Af]     =   fn_forecast_ESN(t, Af, dt_f, Filter);
            tj          =   tj + 1;  
            Aa_extend(:,:,tj)   =   Af;
            t_extend(tj)        =   t;
        end   
        % Store
        ti          =   ti + 1;    
        Aa_KF(:,:,ti:ti+tj-1)	=   repmat(Aa_extend,m,1);
        t_KF(ti:ti+tj-1)      	=   t_extend;
        % Forecast bias in closed loop too
        t_esn  	=   t_bias(end):dt_esn:t_KF(end);% current time to next DA
        Nti     =   length(t_esn);
        [U, ~]  =   advance_ESN('closed', ESN_P, ra, U(:,end), Nti);
        t1_esn  =   t2_esn + 1;
        t2_esn  =   t1_esn + Nti - 2;
        bias(:, t1_esn:t2_esn) 	=   U(:,2:end);
        t_bias(t1_esn:t2_esn)   =   t_esn(2:end);            
    end
end

% ======================================================================= %    
% ======================================================================= % 
