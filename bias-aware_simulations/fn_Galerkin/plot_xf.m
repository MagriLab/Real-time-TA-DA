function plot_xf(Sim)
    if isfield(Sim, 'psi')
        t       =   Sim.t;
        psi     =   Sim.psi;
        P       =   Sim.P;
        u_f     =   P.cos_omjxf * psi(:,1:P.N_m)';   
        p_f     =   - P.sin_omjxf * psi(:,P.N_m+1:2*P.N_m)'; 
        tau_i 	=   round(P.tau / P.dt);
        law     =   P.law;
        beta    =   P.beta;
        kappa   =   P.kappa;
        u_0     =   P.Mean.u_0;
    else
        [p_f, u_f] = deal(zeros(1, length(Sim.ts)));
        for i = 1:length(Sim.ts)
            ti = Sim.ts(i);
            [p_f(i), u_f(i)] = fn_p_from_GH(ti, 0, Sim.Mean, Sim.Geom, ...
                                            Sim.ts, Sim.Gs, Sim.Hs);
        end
        t       =   Sim.ts;
        tau_i   =   round(Sim.Mean.Tau/(t(2)-t(1)));
        law     =   Sim.Mean.Heat_law;
        if contains(law, 'tan')
            beta    =   Sim.Mean.beta;
            kappa   =   Sim.Mean.kappa;
        end
        u_0     =   Sim.Mean.u1;
    end
    u_ftau	=   zeros(size(u_f));   u_ftau(tau_i:end) = u_f(1:end-tau_i+1);
    p_ftau	=   zeros(size(p_f));   p_ftau(tau_i:end) = p_f(1:end-tau_i+1);

    if contains(law, 'G_eqn')
        qdot = Sim.Qs;
    else
        if contains(law, 'tan')
            qdot    =   - (beta^2 * p_ftau) ./ (beta + kappa * u_ftau.^2);
        else
            qdot    =	beta * (sqrt(abs(u_0/3 + u_ftau)) - sqrt(u_0/3)) ; %[kg/s3] = [W/m2]
        end
    end

    figure('Units', 'Normalized', 'OuterPosition', [0 0 0.6 0.5]);
    tiledlayout(2,2)
    nexttile; hold on
    plot(t, p_f)
    xlim([t(1) t(end)])
    xlabel('$t$ [s]'); ylabel('$p''(x_\mathrm{f}, t)$ [Pa]')
    nexttile(3)
    plot(t, p_f)
    xlabel('$t$ [s]'); ylabel('$p''(x_\mathrm{f}, t)$ [Pa]')
    xlim([0.96 0.9998])
    nexttile([2 1]); hold on
    j = find('_' == law);
    if ~isempty(j)
        title([law(1:j-1), ' ',law(j+1:end), ' model'])
    else
        title([law, ' model'])
    end
    plot(t(end), qdot(end))
    plot(t, qdot)    
    xlabel('$t$ [s]'); ylabel('$\dot{q}(t)$ [W/m$^2$]')
end