function [ALL_S, ALL_V, ALL_P] = DA_paper_data(choose_fig)
    % This function creates the datasets required to replicate the
    % figures in Nóvoa, A. & Magri, L. (2021). Real-time thermoacoustic 
    % data assimilation.
    % Input: 
    %   - choose_fig: figure number
    % Output: 
    %   - ALL_S: cell containing the ensemble solution structures. Each 
    %   structure contains the Aa matrix for the true, comparison and 
    %   filtered solutions so as the time vectors and covariance matrix.
    %   - ALL_V: cell containing structures resulting from applying the
    %   obtain_sol_vecs.m to ALL_S
    %   - ALL_P: cell containing the parameters of each solution
    % ===================================================================
    % Define general parameters for the figure
    param  =   DA_paper_params(choose_fig);
    
    % Group the figures according to their type
    modes_figs      =   [5, 7, 9];              % Non-chaos timeseries
    inflation_figs 	=   [11,12];                % QP mics inflation
    chaos_figs   	=   [13,14,16,17,18,19];	% Chaottic figures
    error_figs    	=   [6,10,15];              % QP mics error metrics
    
    % Set the initial and final times of integration as well as the 
    % initial and final times of filtering
    betas	=   [0.2; 0.4; 7.7; 3.6; 7.0];
    t1_vec	=   [0; 900; 900; 900; 900];    % Initial filtering time
    t2_vec	=   [50; 950; 950; 950; 950];   % Final filtering time
    t2A_vec	=   [60; 960; 960; 960; 960];   % Final integration time
    if any(choose_fig == [9,10,13,15,16])
        t2_vec	=   t2_vec + 20; %[70; 970; 970; 970; 970];
        t2A_vec	=   t2_vec + 10; %[80; 980; 980; 980; 980];
    elseif choose_fig == 12
        t2_vec(4)   =   980;	t2A_vec(4)  =   980;
    elseif choose_fig == 11
        t2_vec(4)   =   1000;	t2A_vec(4)  =   1000;
    elseif any(choose_fig == [17,14])
        t2_vec(5)   =   1200;	t2A_vec(5)  =   1500;
    end
    %% ----------------------------------------------------------------- %%
    if any(choose_fig == modes_figs)
        %% NON-CHAOTIC FIGURES
        % Create the true, unfiltered, and filtered solutions for the four 
        % non-chaotic regimes dynamical_modes =   ['FP';'LC';'FL';'QP']; 
        for  vi = 1:4            
            selected_mode  	=   vi;     % Mode to assimilate?
            % Define more parameters
            param.beta   	=   betas(selected_mode);
            param.t_max     =   t2A_vec(selected_mode); %
            param.t_minKF   =   t1_vec(selected_mode);  
            param.t_maxKF   =   t2_vec(selected_mode);  
            param           =   define_parameters(param);
            % REFERENCE SOLUTION
            SOLS                =   [];
            [t_true, Aa_true]   =   create_reference_solution(param);
            SOLS.Aa_true        =   Aa_true; 
            SOLS.t_true         =   t_true; 
            % FALMAN FILTER SOLUTION 
            [t_KF,Aa_KF,Cpp_tr] =   EnSRKF_parfor(param,Aa_true,t_true);
            SOLS.Aa             =	Aa_KF;  
            SOLS.t              =   t_KF;
            SOLS.Cpp_tr         =   Cpp_tr;
            % COMPARISON SOLUTION
            [t_comp, Aa_comp]   =   create_comparison_solution(param,Aa_KF);
            SOLS.Aa_comp        =   Aa_comp;
            SOLS.t_comp         =   t_comp;   
            % OBTAIN SOL VECS
            VECS            =   obtain_sol_vecs(SOLS,param);  
            % STORE FOR LOOP
            ALL_S{vi}   =   SOLS;   
            ALL_V{vi}   =   VECS;
            ALL_P{vi}   =   param;
        end
    elseif any(choose_fig == inflation_figs)
        %% INCREASE, REJECT, INFLATION 
        selected_mode  	=   4;     % QUASIPERIODIC
        % Define pParameters
        param.beta   	=   betas(selected_mode);
        param.t_max     =   t2A_vec(selected_mode); %
        param.t_minKF   =   t1_vec(selected_mode);  
        param.t_maxKF   =   t2_vec(selected_mode);  
        param           =   define_parameters(param);        
        % REFERENCE SOLUTION
        [t_true, Aa_true]   =   create_reference_solution(param);        
        ALL_S   =   {};	ALL_P   =   {};	ALL_V   =   {}; vi = 0;
        if choose_fig == 11
            rho_vec = 1;
        else
            rho_vec = [1 1.02]; 
        end
        for  ri = rho_vec
            param.inflation   =   ri;
            for mi = [10, 50, 100,150]
                vi              =   vi + 1;
                param.m         =   mi;                
                SOLS        	=   [];
                % FALMAN FILTER SOLUTION 
                [t_KF, Aa_KF,~] =   EnSRKF_parfor(param, Aa_true, t_true);
                % KEEP BETA AND TAU ONLY
                SOLS.Aa         =	Aa_KF(:,end-1:end,:);  
                SOLS.t          =   t_KF;                
                % STORE FOR LOOP
                ALL_S{vi}   =   SOLS;
                ALL_P{vi}   =   param;
            end
        end
    elseif any(choose_fig == error_figs)
        %% QUASIPERIODIC ERROR FIGURES
        if choose_fig == 15
            selected_mode   =   5;
        else
            selected_mode  	=   4;     % QUASIPERIODIC
        end
        % Define pParameters
        param.beta   	=   betas(selected_mode);
        param.t_max     =   t2A_vec(selected_mode); %
        param.t_minKF   =   t1_vec(selected_mode);  
        param.t_maxKF   =   t2_vec(selected_mode);  
        param           =   define_parameters(param);
        % REFERENCE SOLUTION
        [t_true, Aa_true]   =   create_reference_solution(param);
        % LEFT HAND SIZE PARAMETER ITERATIONS        
        if choose_fig == 6
            left_vec   =   [0.5, 0.25, 0.1]; % modes std
        elseif choose_fig  == 10
            left_vec   =   [0.1, 0.01, 0.001]; % microphone std
        elseif choose_fig == 15
            left_vec   =   [1.5, 1, 0.5]; % kmeas            
        end
        vi  =   0;
        for  left_i = left_vec
            vi              =   vi + 1;            
            if choose_fig == 6
                param.sig           =   left_i; 
            elseif choose_fig  == 10
                param.mic_sig_frac	=   left_i;
            elseif choose_fig == 15
                param.k_meas   =   left_i / param.dt ;
            end            
            SOLS        	=   [];
            SOLS.Aa_true   	=   Aa_true; 
            SOLS.t_true    	=   t_true; 
            % FALMAN FILTER SOLUTION 
            [t_KF,Aa_KF,Cpp_tr] =   EnSRKF_parfor(param, Aa_true, t_true);
            SOLS.Aa             =	Aa_KF;  
            SOLS.t              =   t_KF;
            SOLS.Cpp_tr         =   Cpp_tr;  
            % OBTAIN SOL VECS
            VECS                =   obtain_sol_vecs(SOLS,param);  
            % STORE FOR LOOP
            ALL_S{vi}   =   SOLS;
            ALL_V{vi}   =   VECS;
            ALL_P{vi}   =   param;
        end
        if choose_fig ~= 15
            % RIGHT HAND SIZE PARAMETER ITERATIONS   
            if choose_fig == 6
                param.sig   =   0.25;
                right_vec   =   [4, 10, 50]; % Ensemble size
            else
                param.mic_sig_frac   =   0.01;
                right_vec   =   [2, 1.5, 1.25]; % kmeas
            end
            for  right_i = right_vec  
                vi              =   vi + 1;
                if choose_fig == 6
                    param.m         =   right_i; 
                else
                    param.k_meas   =   right_i / param.dt ;
                end            
                SOLS        	=   [];
                SOLS.Aa_true   	=   Aa_true; 
                SOLS.t_true    	=   t_true; 
                % FALMAN FILTER SOLUTION 
                [t_KF, Aa_KF, Cpp_tr]       =   EnSRKF_parfor(param, Aa_true, t_true);
                SOLS.Aa         =	Aa_KF;  
                SOLS.t          =   t_KF;
                SOLS.Cpp_tr     =   Cpp_tr;  
                % OBTAIN SOL VECS
                VECS            =   obtain_sol_vecs(SOLS,param);  
                % STORE FOR LOOP
                ALL_S{vi}   =   SOLS;
                ALL_V{vi}   =   VECS;
                ALL_P{vi}   =   param;
            end
        end
    elseif any(choose_fig == chaos_figs)
        %% CHAOTIC FIGURES COMPARING MODES AND MICS ASSIMILATION
        % Create the true, unfiltered, and filtered solutions for the  
        % chaotic assimilation of acoustic Galerkin modes and pressure data
        selected_mode           =   5;     % Chaos mode
        param.beta              =   betas(selected_mode);
        param.t_max             =   t2A_vec(selected_mode); %
        param.t_minKF           =   t1_vec(selected_mode);  
        param.t_maxKF           =   t2_vec(selected_mode);
        param.assimilate_modes  =   false;
        param                   =   define_parameters(param);
        % REFERENCE SOLUTION
        SOLS                    =   [];
        [t_true, Aa_true]       =   create_reference_solution(param);
        SOLS.Aa_true            =   Aa_true; 
        SOLS.t_true             =   t_true;
        ai                      =   0;     
        if choose_fig ~= 14
            for a = [true, false]
                ai  =    ai + 1;
                % Define parameters
                param.assimilate_modes  =   a;
                param                   =   define_parameters(param);
                % FALMAN FILTER SOLUTION 
                [t_KF,Aa_KF,Cpp_tr]	=   EnSRKF_parfor(param,Aa_true,t_true);
                SOLS.Aa            	=	Aa_KF;  
                SOLS.t             	=   t_KF;
                SOLS.Cpp_tr        	=   Cpp_tr;
                % COMPARISON SOLUTION
                [t_comp, Aa_comp]   =   create_comparison_solution(param,Aa_KF);
                SOLS.Aa_comp        =   Aa_comp;
                SOLS.t_comp         =   t_comp;   
                % OBTAIN SOL VECS
                VECS                =   obtain_sol_vecs(SOLS,param);
                % STORE FOR LOOP
                if choose_fig == 17
                    VECS  	=   rmfield(VECS,{'eta_true','eta','eta_comp',...
                                        'v_xf','v_xf_true','v_xf_comp',...
                                        'etad_true','etad','etad_comp'});
                    SOLS 	=   rmfield(SOLS, {'Aa','Aa_comp','Cpp_tr'});
                end                
                ALL_S{ai}   =   SOLS;
                ALL_V{ai} 	=   VECS;
                ALL_P{ai} 	=   param;  
            end
        else
%             VECS     	=   obtain_sol_vecs(SOLS,param);
            ALL_S{1}    =   SOLS;
            ALL_V{1}    =   [];
            ALL_P{1}	=   param;
        end
        
    end

end
