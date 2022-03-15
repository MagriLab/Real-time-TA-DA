function fn_check_esn(Sim_t, Sim_f, bias, t_bias, ESN_P)
    N_m         =   Sim_f.P.N_m;
    try
        t1          =   1/Sim_t.P.dt;
    catch
        t1          =   1/(Sim_t.t_mic(2)-Sim_t.t_mic(1)) ;
    end
    t_interval  =   t1:t1+500;
    t_t         =   Sim_t.t_mic(t_interval);
    p_ht        =   Sim_t.p_mic(:,t_interval);
    sins        =   sin(Sim_f.P.omega_j./Sim_f.P.Mean.c_0 * Sim_f.P.x_mic)'; 
    p_hf        =   - sins * Sim_f.psi(:,N_m+1:2*N_m)'; 
    t_f         =   Sim_f.t;
    % ------------------------- Create washout -------------------------- %
    N_wash      =   ESN_P.N_washout;
    dt_esn      =   ESN_P.dt;
    [U_wash, t_wash]	=   create_washout(dt_esn,N_wash,t_f,p_hf,t_t,p_ht);

    % ------------------- Initialise ESN in open loop ------------------- %
    [U_open, ra] =   advance_ESN('open', ESN_P, U_wash);

    U0   =   U_wash(end,:);
    % ------------------- Forecast ESN in closed loop ------------------- %
    t_esn   =   t_t(end):dt_esn:t_t(end)+0.5;% From current time to next DA
    Nti_esn =   length(t_esn);
    [U, ~]  =   advance_ESN('closed', ESN_P, ra, U0, Nti_esn);

    % -------------------- Plot bias and ESN output --------------------- %
    figure; hold on
    l1=plot(t_bias, bias(1,:),'LineWidth',3); l1.Color(4)=0.5;
    l=plot(t_esn, U(1,:),'LineWidth',2); c=get(l, 'color');
    plot(t_wash, U_open(:,1), 'o--','color',c)
    plot(t_wash, U_wash(:,1), 'x--','color','k')
    xlim([t_wash(1)-0.01, t_wash(end)+0.05])
    legend('actual bias','closed loop','open loop','washout')
end