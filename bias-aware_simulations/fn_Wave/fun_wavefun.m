function [gt,ht,gt_Tu,Qt] = fun_wavefun(varargin)%t,ts,Xis,Gs,Hs,j,J,Mean,Geom)
    % Function that computes g(t) and h(t) via Eq (3.9)
    if nargin == 9
        [t,ts,Xis,Gs,Hs,j,J,Mean,Geom]=deal(varargin{:});
    elseif nargin == 7
        [mi,t,Sim,j,J,Mean,Geom]=deal(varargin{:});        
        Gs = Sim.Gs(mi,:);
        Hs = Sim.Hs(mi,:);
        ts = Sim.ts;
        Xis = Sim.Xis;
    else
        error('wrong number of inputs to fun_wavefun.m')
    end
    % Retrieve variables:
    c1 = Mean.c1;
    u1 = Mean.u1;
    M1 = Mean.M1;
    rho1 = Mean.rho1;
    R_in    =   Mean.R_in;
    R_out   =   Mean.R_out;
    % Network paramters:
    Tu = Mean.Tu;
    Td = Mean.Td;
    InvX  = Mean.IX;
    InvXY = Mean.IXY;
    % Flame parameters:
    Abar = Mean.Abar;
    Qbar = Mean.Qbar;
    Tau  = Mean.Tau;
    % Geometric parameters
    A1 = Geom.A1;
    r  = Geom.r;
    % ============== Compute unsteady heat release ===================== %

    if ~isfield(Mean, 'Heat_law') || strcmp(Mean.Heat_law,'G_eqn')    
%         if Qbar == 0
%             Qt = 0;
%         else
            % Retrieve wave amplitudes at previous states:
            gt_Tau    = fwavfun_hist(t-Tau,Gs,ts,j,J);
            gt_Tu_Tau = fwavfun_hist(t-Tu-Tau,Gs,ts,j,J);
            gt_Tu     = fwavfun_hist(t-Tu,Gs,ts,j,J);
            ht_Td     = fwavfun_hist(t-Td,Hs,ts,j,J);
            % Flame position at previous states:
            xit_Tau = fn_flame_front_hist(t-Tau,Xis,ts,j,J,Mean);
            % Velocity U(t - Tau)
            Ut_Tau = fn_flame_holder_vel(t,gt_Tau,gt_Tu_Tau,Mean,Geom);
            % dxi/dr(t - Tau)
            dxi_dr = fn_dxi_dr(xit_Tau,Ut_Tau,Mean,Geom);
            % Computes flame Area: A(t - Tau) Eq. (2.9)
            A_t_Tau = trapz(r,2*pi*r.*sqrt(1 + dxi_dr.^2));
            % Compute Heat release Q(t). Eq. (3.1)
            Qt = Qbar*A_t_Tau/Abar;
%         end
        Fq = [0;(Qt - Qbar)/(A1*c1)];
    else 
        beta   	=   Mean.beta;
        kappa 	=   Mean.kappa;
        
        st  =   max(1,j-J);
        g   =   @(t) fv_interp(t,ts(st:j-1),Gs(st:j-1));
        f   =   @(t) R_in * g(t - Tu);
        F   =   f(t-Tau); % at flame loc thi si F = f(t-Tau)
        G   =   g(t-Tau);
        

        u_ftau  =   1/(rho1*c1)*(F - G);
        p_ftau  =   (F + G);

        
        if t < Tau
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
        
%         Qt	=   ;
        
% 
        end


        gt_Tu	=   fwavfun_hist(t-Tu,Gs,ts,j,J);
        ht_Td	=   fwavfun_hist(t-Td,Hs,ts,j,J);

    end
    
    % Computes g(t) and h(t) via Eq. (3.9), where
    % X*w(t) = Y*w(t-Tau) + Fr(t), w = [g(t);h(t)]

%     Fq = [0;(Qt - Qbar)/(A1*c1)];
    
    w = InvXY * [gt_Tu;ht_Td] + InvX * Fq;
    % Output:
    gt = w(1);
    ht = w(2);
end

% Function definitions:
function v = fv_interp(t,ts,Vs)
    % Interpolate
    if t < 0
        v = 0;
    else
        v =  interp1(ts,Vs,t);
    end
end