function [t_KF, Aa_KF, Obs, t_bias, bias] =   EnSRKF_Wave(Truth, Filter, ESN_P)

    % Retrieve true variables
    t_true  =   Truth.t_mic;
    p_true  =   Truth.p_mic;
    % Retrieve filter variables
    t_max   =   Filter.t_max;
    t_stop  =   Filter.t_stop;
    m       =   Filter.m;
    N       =   Filter.N;
    N_m     =   Filter.N_m;
    N_mic   =   Filter.N_mic;
    dt      =   Filter.dt;
    
    
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
    % Select 10 points of the history
    N_h     =   10;
    x_mic   =   Filter.x_mic;
    Mean    =   Filter.Mean;
    pA      =   zeros(m,N_mic);
    % ====================== Loop over observations ===================== % 
    while true
        % ----------------- Create state vector/matrix ------------------ %
        % Af = [Gs(t),...,Gs(t-J),Hs(t),...,Hs(t-J),params]'
        for j = 1:m
            tau     =   Af(j,end);
            J       =   ceil(max([1,[tau, tau+Mean.Tu, tau+Mean.Td]/dt]));
            G_hist  =   squeeze(Aa_KF(j,1,end-J+1:end));
            H_hist  =   squeeze(Aa_KF(j,N_h+1,end-J+1:end));
            t_hist  =   t_KF(end-J+1:end);
            % Fit the history to N_h points
            t_Nh    =   linspace(t_hist(1),t_hist(end),N_h);
            G_Nh    =   interp1(t_hist, G_hist, t_Nh);
            H_Nh    =   interp1(t_hist, H_hist, t_Nh);
            % Store into Af
            Af(j,1:2*N_h)   =   [G_Nh, H_Nh];            
            % Compute biased pressure state
            [~,pA(mi,:)]	=   fn_p_from_GH(t,x_mic,Mean,t_hist,G_hist,H_hist);
        end
        % ----------- Augment forecast with unbiased pressure ----------- %
        pA_m  	=   mean(pA,2);
        pA   	=   (pA - pA_m) + (pA_m + U); % Deviations + updated mean
        Af_aug	=   vertcat(Af',pA); % [m x (N + N_mic)]
        % ------------ Perform batch initialisation/analysis ------------ %%%%%% Continue here%
        [Aa, py, Filter] = fn_analysis_wave(ti_mic, Af_aug, p_mic, Filter);
        % Update solution at time t
        Aa_KF(:,:,ti)	=   Aa;
        Obs(:,ti_mic)	=   py;
        t_Obs(ti_mic)   =   t_mic;
        % get next time for analysis and break if no more data is available
        try 
            ti_mic  =   ti_mic + 1;  
            t_mic   =   t_true(ti_mic);
            p_mic   =   p_true(ti_mic,:); 
            Af      =   Aa;
        catch
            break
        end  
        % ================ Loop over forecast time steps ================ % 
        stop_f = false;
        while ~stop_f
            if	t + dt < t_mic
                dt_f    =   dt;
            else
                dt_f    =   t + dt - t_mic;
                stop_f	=   true;
            end
            % ------------ Forecast members individually ------------ % 
            [t, Af]     =   fn_forecast_ESN(t, Af, dt_f, Filter);
            ti          =   ti + 1;        
            % Store
            Aa_KF(:,:,ti)	=   Af;
            t_KF(ti)        =   t;
        end       
        % ================== Loop over ESN time steps =================== % 
%         stop_esn = true;
%         while ~stop_esn
            ESN_start_i     =   50;
            if ti_mic >= ESN_start_i
                if ti_mic == ESN_start_i
                    % ----------- Initialise ESN in open loop ----------- %
                    p_h_t 	=   Obs;
                    p_h_KF	=   - sin_omj * ...
                                squeeze(mean(Aa_KF(:,N_m+1:2*N_m,:),1));
                    U_wash	=   create_washout(dt_esn, N_washout, ...
                                t_KF, p_h_KF, t_Obs, p_h_t);
                    [U, ra] =   advance_ESN('open', ESN_P, U_wash);
                else
                    % -------------- Update bias with data -------------- %
%                     U   =   py' - pA_m;
                end
                % ------------- Forecast ESN in closed loop ------------- %
                t_esn   =   t:dt_esn:t_mic;% From current time to next DA
                Nti_esn =   length(t_esn);
                [U, ra] =   advance_ESN('closed', ESN_P, ra, U, Nti_esn);
                % Store
                t1_esn  =   length(t_bias) + 1;
                t2_esn  =   t1_esn + Nti_esn -1;
                
                bias(:, t1_esn:t2_esn) 	=   U;
                t_bias(t1_esn:t2_esn)	=   t_esn;
                
                % Keep only last element, i.e. U(t_mic)
                U   =   U(:,end);
            end
%             t_esn   =   te;
%         end
            
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

function U_washout = create_washout(dt_ens, N_washout, t_KF, p_KF, t_true, p_true)
    % Define washout time
    t_washout   =   t_true(1):dt_ens:t_true(end);
    % Evaluate pressure at washout timesteps
    p_true      =   interp1(t_true, p_true', t_washout);
    p_KF        =   interp1(t_KF, p_KF', t_washout);
    % Bias
    U_washout   =   p_true - p_KF;
    U_washout   =   U_washout(end-N_washout+1:end,:);
end