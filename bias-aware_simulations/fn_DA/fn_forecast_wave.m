function [Af, t] = fn_forecast_wave(t, Af, dt_f, Filter)
	% Retrieve variables
    N_h     =   Filter.N_h;
    Tau     =   Mean.Tau;
    J       =   ceil(max([1,[Tau, Tau+Mean.Tu, Tau+Mean.Td]/dt]));
    % =================== Loop over ensemble members ==================== %
    for j = 1:size(Af,1)
        Gs      =   Af(j, 1:N_h);
        Hs      =   Af(j, N_h+1:2*N_h);
        tau     =   Af(j, end);
        ts      =   linspace(t-tau, t, N_h);
        % Define g and f functions at flame location
        g           =   @(t) fn_hist(t,ts,Gs,J);
        f           =   @(t) R_in * g(t - Tu);
        gt_Tau      = 	g(t-Tau);
        ft_Tau      = 	f(t-Tau);
        if ~isfield(Mean, 'Heat_law') || strcmp(Mean.Heat_law,'G_eqn')    
            % ========================================================== %
            A1      =   Geom.A1;
            u1      =   Mean.u1;
            c1      =   Mean.c1;
            rho1    =   Mean.rho1;
            R_in    =   Mean.R_in;
            Qbar    =   Mean.Qbar;
            if Qbar == 0
                Qt = 0;
            else
                % ---------------------- G equation --------------------- %
                % Compute the velocity U(t - Tau) = 1/(rho1 * c1)  * 
                % (f(t-Tau) - g(t-Tau)) w/ f(t-Tau) = Rin * g(t-Tau-tau_u)
                xit_Tau     =   fn_hist(t-Tau,Xis,ts,J);       
                u_ftau      =   u1 + 1/(rho1*c1) * (ft_Tau - gt_Tau);                 
                % Compute flame aerea as a function of dxi/dr(t - Tau)
                dxi_dr      =   fn_dxi_dr(xit_Tau,u_ftau,Mean,Geom);
                A_t_Tau     =   trapz(r,2*pi*r.*sqrt(1 + dxi_dr.^2));
                % Compute Heat release Q(t). Eq. (3.1)
                Qt          =   Qbar * A_t_Tau/Abar;
            end
            Fq = [0;(Qt - Qbar)/(A1*c1)];
        else % ========================================================== %
            beta   	=   Mean.beta;
            kappa 	=   Mean.kappa;



            u_ftau  =   u1 + 1/(rho1*c1)*(ft_Tau - gt_Tau);
            p_ftau  =   p1 + (F + G);


            if t < tau
                u_ftau = .1;    p_ftau = -.1;
            end
            % Switch x axes for the Q laws
            L       =   Geom.Lb + Geom.Lu;
            xf      =   Geom.Lu; 
            if contains(Mean.Heat_law, 'tan')
            % --------------------------- atan ------------------------------ %
                Qt	=   - (beta^2 * p_ftau) / (beta + kappa * u_ftau^2) * ...
                        cos(pi*xf/L) / sin(pi*xf/L); % [Pa/s]. Note [beta] = 1/s [kappa] = s/m2

        %         Qt	=   beta * sqrt(beta/kappa) * ...
        %                  	atan(sqrt(kappa/beta) * u_ftau); %[m/s2] 

                Fq 	=   [0; Qt / (pi*c1/L) ];
            elseif contains(Mean.Heat_law, 'sqrt')
            % --------------------------- sqrt ------------------------------ %
                Qt	=	beta * (sqrt(abs(u1/3 + u_ftau)) - sqrt(u1/3)) ; %[kg/s3] = [W/m2]
                Qt  =   Qt * A1; %[W]
                Fq 	=   [0;(Qt - Qbar)/(A1*c1)]; %[Pa]

            elseif contains(Mean.Heat_law, 'cub')
            % --------------------------- cubic ----------------------------- %
            % TBD
            end
            % ========================================================== %
        end
        
        gt_Tu       = 	fn_hist(t-Tu,Gs,ts,J);
        ht_Td       =	fn_hist(t-Td,Hs,ts,J);

        % Compute g(t) and h(t) via Eq. (3.9), where w = [g(t);h(t)] and
        % X*w(t) = Y*w(t-Tau) + Fr(t)

        w = InvXY * [gt_Tu;ht_Td] + InvX * Fq;
        % Output:
        gt = w(1);
        ht = w(2);
        
        
        
        
        % ----------------- Create state vector/matrix ------------------ %
        % Af = [Gs(t),...,Gs(t-J),Hs(t),...,Hs(t-J),params]'
        G_hist  =   squeeze(Aa_KF(j,1,end-J+1:end));
        H_hist  =   squeeze(Aa_KF(j,N_h+1,end-J+1:end));
        t_hist  =   t_KF(end-J+1:end);
        % Fit the history to N_h points
        t_Nh    =   linspace(t_hist(1),t_hist(end),N_h);
        G_Nh    =   interp1(t_hist, G_hist, t_Nh);
        H_Nh    =   interp1(t_hist, H_hist, t_Nh);
        % Store into Af
        Af(j,1:2*N_h)   =   [G_Nh, H_Nh];            
    end
end

