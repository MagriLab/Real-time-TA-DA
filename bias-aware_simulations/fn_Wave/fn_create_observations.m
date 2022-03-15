function T = fn_create_observations(model)

    if contains(model, 'wave','IgnoreCase',true); fprintf('Wave approach:')  
        % =============================================================== %        
        if contains(model, 'low'); %fprintf(' low order model ')
        else; %fprintf(' high order model ')
        end
        % Setup the simulation
        Sim     =   setup_sim(model);
        % Time march the simulation. This is an origunal code.
        Geom   	=   Sim.Geom;
        Mean	=   Sim.Mean;
        t_end   =   Sim.Measurement.Time;
        [ts,Gs,Hs,Xis,Qs] = fn_time_marching(t_end,Sim.dt,Geom,Mean); 
        % Compute pressure at microphone location at the measurement frequency
        [t_mic, p_mic, u_mic] = GH_to_pu(Sim,ts,Gs,Hs);

        % Store true data
        T.p_mic         =   p_mic;
        T.u_mic         =   u_mic;
        T.t_mic         =   t_mic;
        T.Gs            =   Gs;
        T.Hs            =   Hs;
        T.ts            =   ts;
        T.Xis           =   Xis;
        T.Qs            =   Qs;
        T.Geom          =   Geom;
        T.Mean          =   Mean;
        T.Measurement	=   Sim.Measurement;
        T.Forcing       =   Sim.Forcing;
        % =============================================================== %
    elseif contains(model, 'Galerkin','IgnoreCase',true); %fprintf('Galerkin approach ')  
        % =============================================================== %        
        Sim =   setup_sim(model);
        P   =   setup_P_dim(Sim);
        
        [t_vec, psi_vec]	=   ode45(@(t,y) gov_eqns_dim(t,y,P,P.law), ...
                                  (P.t_min:P.dt:P.t_max), P.IC);
        
%         ALPHAj = psi_vec(:,P.N_m+1 : 2*P.N_m);
%         p_f     =   - P.sin_omjxf * ALPHAj';   % Pressure at the flame at t
%         
%         figure;
%         plot(t_vec, p_f)


%         % Integrate dim equations
%         t_ode       =   t_min:P.dt:t_max;
%         psi  	=   P.IC; t = t_ode(1);   
%         psi_vec	=   zeros([length(psi),length(t_ode)]);
%         t_vec   =   zeros([1,length(t_ode)]);
%         for ti = 1:length(t_ode) 
%             psi_vec(:,ti)	=   psi;
%             t_vec(ti)  	=   t;
%             fun         =   @(t,psi) gov_eqns_dim(t,psi,P,P.law);
%             [t, psi]    =   RK4_solver(t,psi,P.dt,fun);
%         end
        
        % Compute pressure and velocity at the microphones' locations
        [t_mic, p_mic, u_mic] = psi_to_pu(t_vec,psi_vec,P);
        
        if size(psi_vec,1) < size(psi_vec,2)
            psi_vec = psi_vec';
        end
        
        % Output
        T.t     =   t_vec';
        T.psi	=   psi_vec;
        T.P     =	orderfields(P);        
        T.p_mic	=   p_mic';
        T.u_mic	=   u_mic';
        T.t_mic =   t_mic;
        % =============================================================== %
    else
        error('Type of solver not defined')
    end
    T	=	orderfields(T);
end

% ======================================================================= %




% ----------------------------------------------------------------------- %
