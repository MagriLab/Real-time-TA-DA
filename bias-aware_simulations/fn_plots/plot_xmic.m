function plot_xmic(t, psi, P, x_mic)

    u_f     =   P.cos_omjxf * psi(:,1:P.N_m)';   
    p_f     =   - P.sin_omjxf * psi(:,P.N_m+1:2*P.N_m)';  
    
    tau_i 	=   P.tau / P.dt; 
    u_ftau	=   zeros(size(u_f));   u_ftau(tau_i:end) = u_f(1:end-tau_i+1);
    p_ftau	=   zeros(size(p_f));   p_ftau(tau_i:end) = p_f(1:end-tau_i+1);

    if contains(P.law, 'tan')
        qdot    =   - (P.beta^2 * p_ftau) ./ (P.beta + P.kappa * u_ftau.^2);
    else
        qdot    =	P.beta * (sqrt(abs(P.u_0/3 + u_ftau)) - sqrt(P.u_0/3)) ; %[kg/s3] = [W/m2]
    end


    figure('Units', 'Normalized', 'OuterPosition', [0 0 0.8 0.7]);
    tiledlayout(length(x_mic)/2, 6, 'Padding', 'compact', 'TileSpacing',  'compact')
    for i = 1:length(x_mic)
        p_f	=   - sin(P.omega_j'./P.c_0 * x_mic(i)) * psi(:,P.N_m+1:2*P.N_m)'; 
        nexttile; hold on
        plot(t, p_f)
        xlabel('$t$ [s]'); ylabel('$p''(x_\mathrm{f}, t)$ [Pa]')
        nexttile;
        plot(t, p_f)
        xlabel('$t$ [s]');
        legend(['$p''(x=$',num2str(x_mic(i)),' m $, t)$'])
        xlim([0.96 01])
        if i ==length(x_mic)
            nexttile([1,3]); hold on
            plot(t, qdot);plot(t, qdot)
            xlabel('$t$ [s]');
            ylabel('$\dot{q}''(t)$ [W/m$^2$]')
            nexttile([1,3]); hold on
            plot(t, qdot);plot(t, qdot)
            xlim([0.96 1])
            xlabel('$t$ [s]');
            Z = 1000;
            ylabel('$\dot{q}''(t)$ [W/m$^2$]')
        end
    end

end