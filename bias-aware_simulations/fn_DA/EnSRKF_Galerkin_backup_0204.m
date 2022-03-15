function [t_KF, Aa_KF, Obs, t_bias, bias] =   EnSRKF_Galerkin_backup_0204(Truth, Filter, ESN_P)

%  THIS VERSION APPLIES OK THE ESN WITH ASSIMILATION, ESN AND FORECAST dts
%  MULTIPLES OF EACH OTHER. NO DATA ASSIMILATION AFTER AN INITIAL
%^%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




    % Retrieve true variables
    t_true  =   Truth.t_mic;
    p_true  =   Truth.p_mic;
    % Retrieve filter variables
    t_max   =   Filter.t_max;
    t_stop  =   Filter.t_stop;
    t_start =   Filter.t_start;
    m       =   Filter.m;
    N       =   Filter.N;
    N_m     =   Filter.N_m;
    N_mic   =   Filter.N_mic;
    dt      =   Filter.dt;
    dt_f    =   dt;
    
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
    U       =   zeros(N_mic, 1);
    Obs     =   zeros(N_mic, 1);
    t_Obs   =   t_mic;
    bias    =   zeros(N_mic, 1);
    t_bias  =   [];
    
    % ====================== Loop over observations ===================== % 
    while true
        % -------------- Compute forecast pressure states --------------- %
        pA    	=   - sin_omj * Af(:,N_m+1:2*N_m)';
        pA_m  	=   mean(pA,2);
        pA   	=   (pA - pA_m) + (pA_m + U(:,end)); % Deviations + updated mean
        Af_aug	=   vertcat(Af',pA); % [m x (N + N_mic)]
%         % ------------ Perform batch initialisation/analysis ------------ %
        [Aa, py, Filter] = fn_analysis_ESN(ti_mic, Af_aug, p_mic, Filter);
        % Update solution at time t
        if ti_mic <= 10
            Aa_KF(:,:,ti)	=   Aa;
            Af      =   Aa;
        end
            
        Obs(:,ti_mic)	=   py;
        t_Obs(ti_mic)   =   t_mic;
        % get next time for analysis and break if no more data is available
        try 
            ti_mic  =   ti_mic + 1;  
            t_mic   =   t_true(ti_mic);
            p_mic   =   p_true(ti_mic,:); 
%             Af      =   Aa;
        catch
            break
        end         
        % ================== Loop over ESN time steps =================== % 
%         stop_esn = true;
%         while ~stop_esn
            ESN_start_time     =   0.5;
            if (t_mic - t_start) >= ESN_start_time
                if sum(U) == 0
                    % ----------- Initialise ESN in open loop ----------- %
                    p_ht 	=   Obs;
                    p_hf	=   - sin_omj * squeeze(mean(Aa_KF(:,N_m+1:2*N_m,:),1));
                    [U_wash, t_w]   =   create_washout(dt_esn, 2*N_washout, ...
                                            t_KF, p_hf, t_Obs, p_ht);
                    [U_open, ra]	=   advance_ESN('open', ESN_P, U_wash);
                    U   =   U_open(end,:)';
                    bias(:,1)   =   U;
                    t_bias(1)   =   t-dt; 
% %                 %%%%%%%%%%%% PLOT DEBUG
%                     figure; hold on
%                     plot(t_w, U_open(:,1), 'b*--')
%                     l=plot(t_w, U_wash(:,1), 'k-','linewidth', 3); l.Color(4) = 0.2;
% %                 %%%%%%%%%%%% PLOT DEBUG


                    % -------------- Update bias with data -------------- %
                end
                pA_m	=   - sin_omj * mean(Aa_KF(:,N_m+1:2*N_m,end),1)';
                if  length(t_bias) > 2
                    U       =   py' - pA_m;
                end

                % ------------- Forecast ESN in closed loop ------------- %
                t_esn   =   t_bias(end):dt_esn:t_mic;% From current time to next DA
                Nti_esn =   length(t_esn);
                [U, ra] =   advance_ESN('closed', ESN_P, ra, U, Nti_esn);
                
%                 %%%%%%%%%%%% PLOT DEBUG
%                 p_ht 	=   Obs;
%                 p_hf	=   - sin_omj * squeeze(mean(Aa_KF(:,N_m+1:2*N_m,:),1));
%                 [U_wash, t_w]   =   create_washout(dt_esn, N_washout/2, ...
%                                         t_KF, p_hf, t_Obs, p_ht);
%                 l=plot(t_w, U_wash(:,1), 'k-','linewidth', 3); l.Color(4) = 0.2;
%                 plot(t_esn, U(1,:), 'r*-')
%                 plot(t_Obs(end), py(1) - pA_m(1), 'kx','MarkerSize',12)
%                 %%%%%%%%%%%% PLOT DEBUG
                
                % Store
                t1_esn  =   length(t_bias) + 1;
                t2_esn  =   t1_esn + (Nti_esn - 1) - 1;
                
                bias(:, t1_esn:t2_esn) 	=   U(:,2:end);
                t_bias(t1_esn:t2_esn)	=   t_esn(2:end);
                
                % Keep only last element, i.e. U(t_mic)
%                 U   =   U(:,end);
            end
%             t_esn   =   te;
%         end
        % ================ Loop over forecast time steps ================ % 
        stop_f = false;
        while ~stop_f % ========================== LOOP OVER FORECAST STEPS
%             if	(t + dt - t_mic) < 1e-8
%                 if dt_f ~= dt
%                     dt_f    =   dt - (t - t_KF(ti-1));
%                 else
%                     dt_f    =   dt;
%                 end
%             else
%                 dt_f    =   t + dt - t_mic;
%             end
            dt_f    =   dt;
            if	t + dt >= t_mic
                stop_f	=   true;
            end
            % ------------ Forecast members individually ------------ % 
            [t, Af]     =   fn_forecast_ESN(t, Af, dt_f, Filter);
            ti          =   ti + 1;        
            % Store
            Aa_KF(:,:,ti)	=   Af;
            t_KF(ti)        =   t;
%             if (t - t_mic) < 1e-8
%                 t = t_mic;
%                 stop_f	=   true;
%             end
        end
            
    end

    % ============== Integrate further as mean if required ============== %    
    if t_stop < t_max
        Af  =   mean(Af,1)';
        tj  =   0;
        Aa_extend	=   zeros(1, N, 1);
        t_extend    =   zeros();
        while t < t_max
            [t, Af] = fn_forecast_ESN(t, Af, dt_f, Filter);
            tj          =   tj + 1;  
            Aa_extend(:,:,tj)   =   Af;
            t_extend(tj)        =   t;
        end
        ti          =   ti + 1;        
        % Store
        Aa_KF(:,:,ti:ti+tj-1)	=   repmat(Aa_extend,m,1);
        t_KF(ti:ti+tj-1)      	=   t_extend;
    end

    % ======================= Store output vectors ====================== %    
%     Aa_KF   =   Aa;
%     t_KF    =   TX;
end

% ======================================================================= %    
% ======================================================================= % 
