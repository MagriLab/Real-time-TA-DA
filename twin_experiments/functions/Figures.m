function OUT = Figures(S,V,P,fig_type,folder)
    %  plot_settings
    if ~isempty(folder)
        suffix  = ['_micsig',num2str(P.mic_sig_frac),'_kmeas',...
                    num2str(P.k_meas),'_IC',num2str(P.IC_param_fact)];
    end
    OUT =[];
    if contains(fig_type,'timeseries')% ------------------------------------------
        box on;    hold on        
        if  isfield(V, 'p_xf_comp')
            plot(S.t_true,V.p_xf_true,'linewidth',1.5,'DisplayName','True sol.')
            plot(S.t_comp,V.p_xf_comp,'linewidth',0.5,'DisplayName','Unfiltered sol.')
            plot(S.t,V.p_xf,'linewidth',1,'DisplayName','Filtered sol.')
        else
            if isfield(V, 'p_xf_true')
                plot(S.t_true,V.p_xf_true,'g','linewidth',2.5,'DisplayName','True sol.')
            end
            plot(S.t,V.p_xf, 'k--','linewidth',1,'DisplayName','Filtered sol.')
        end
        for i = 1:P.k_meas:(P.t_maxKF-P.t_minKF)/P.dt+1
            analysis_step_t(i)  =   S.t(i);
            analysis_steps_p(i) =   V.p_xf(i);
        end
        plot(analysis_step_t,analysis_steps_p,'ko','MarkerSize',5,'MarkerFaceColor', 'r', 'DisplayName','Analysis step')
        hold off
        xlim([min(S.t), max(S.t)])
        leg = legend;
        xlabel('$t$');     ylabel('$p''_\mathrm{f}$');
        set(leg,'Orientation','horizontal','Location','northoutside','EdgeColor',[1 1 1])
        if ~isempty(folder)
            savefig(fig1,[folder,'timeseries',suffix,'.fig']);
        end
    end  
    
    if contains(fig_type,'zoom')%--------------------------------------------------
        box on;    hold on        
        if  isfield(V, 'p_xf_comp')
            plot(S.t_true,V.p_xf_true,'linewidth',1.5,'DisplayName','True sol.')
            plot(S.t_comp,V.p_xf_comp,'linewidth',0.5,'DisplayName','Unfiltered sol.')
            plot(S.t,V.p_xf,'linewidth',1,'DisplayName','Filtered sol.')
        else
            if isfield(V, 'p_xf_true')
                plot(S.t_true,V.p_xf_true,'g','linewidth',2.5,'DisplayName','True sol.')
            end
            plot(S.t,V.p_xf, 'k--','linewidth',1,'DisplayName','Filtered sol.')
        end
        for i = 1:P.k_meas:(P.t_maxKF-P.t_minKF)/P.dt+1
            analysis_step_t(i)  =   S.t(i);
            analysis_steps_p(i) =   V.p_xf(i);
        end
        plot(analysis_step_t,analysis_steps_p,'ko','MarkerSize',8,'MarkerFaceColor', 'r', 'DisplayName','Analysis step')
        xlim([max(S.t)-14, max(S.t)-4])
        set(gca, 'FontSize', 24)
        leg = legend;
        set(leg,'Orientation','horizontal','Location','northoutside','EdgeColor',[1 1 1])
        if ~isempty(folder)
            savefig(fig1,[folder,'timeseries',suffix,'.fig']);
        end
    end  
    if contains(fig_type,'params')%-----------------------------------------------

        beta_m = squeeze(S.Aa(:,end-1,:))'/P.beta;
        tau_m = squeeze(S.Aa(:,end,:))'/P.tau;
        
        if  contains(fig_type,'params_2')
            hold on;
            ylabel('$\beta (t) / \beta^\mathrm{true}$');
            dispname    =   ['$m=$', num2str(P.m)]; 
            dispname2   =   ''; 
        else
            dispname    =   '$\beta (t) / \beta^\mathrm{true}$';
            dispname2	=   'std($\beta (t) / \beta^\mathrm{true}$)';
        end

        mean_   =      mean(beta_m)';
        std_p   =      mean_+std(beta_m)';
        std_m   =      mean_-std(beta_m)';
        l   =   plot(S.t(:),mean_,'LineWidth',1.4,'DisplayName', dispname); 
        c   =   get(l,'Color');
        patch([S.t(:); flipud(S.t(:))],[std_m; flipud(std_p)],c,...
            'FaceAlpha',0.2, 'EdgeColor','none','DisplayName',dispname2); 
          
        if  contains(fig_type,'params_2')
            dispname    =   ['$m=$', num2str(P.m)]; 
            dispname2   =   ''; 
            line_style 	=   '--';
            loc         =   'eastoutside';
            col         =   1;
            
            legend('Location',loc, 'NumColumns', col); legend('boxoff')
            nexttile(str2double(fig_type(end))); hold on
            ylabel('$\tau (t) / \tau^\mathrm{true}$');
        else
            c = [0.85 0.325 0.098];
            dispname   =   '$\tau (t) / \tau^\mathrm{true}$';
            line_style  =   '-';
            loc         =   'southeast';
            col         =   2;
        end
          
        mean_   =      mean(tau_m)';
        std_p   =      mean_+std(tau_m)';
        std_m   =      mean_-std(tau_m)';

        hold on;
        plot(S.t(:),mean_,line_style,'color',c, 'LineWidth',1.4,'DisplayName', dispname); 
        hold on;
        patch([S.t(:); flipud(S.t(:))],[std_m; flipud(std_p)],c, 'FaceAlpha',0.2, ...
              'EdgeColor','none','DisplayName',dispname2); 

        legend('Location',loc, 'NumColumns', col)
        xlabel('$t$');
        legend('boxoff')

    end
    if contains(fig_type,'PSD')%-----------------------------------------------
        % Power spetral density comparison true and filtered
        %True solution PSD
        psd_t_d =   0.5 * (P.t_maxKF - P.t_minKF)/ P.dt ;
        ax1 = gca;
        [f,PSD_]    =   fun_PSD(P.dt,V.p_xf_true(P.t_maxKF / P.dt - psd_t_d:P.t_maxKF / P.dt));
        OUT.true_PSD    =   PSD_;
        OUT.true_f    =   f;
        plot(f, PSD_,'LineWidth',2.5,'DisplayName', 'True sol.');
        ylabel('PSD'); hold on;
        
        if  isfield(V, 'p_xf_comp')
            [f,PSD_] =   fun_PSD(P.dt,V.p_xf_comp(end-psd_t_d:end));
            OUT.filter_PSD      =   PSD_;
            OUT.filter_f        =   f;
            plot(f, PSD_, 'LineWidth',1, 'DisplayName', 'Unfiltered sol.');
        end
        
        
        [f,PSD_] =   fun_PSD(P.dt,V.p_xf((P.t_maxKF - P.t_minKF)/ P.dt - psd_t_d:(P.t_maxKF - P.t_minKF)/ P.dt));
        OUT.filter_PSD	=   PSD_;
        OUT.filter_f  	=   f;
        plot(f, PSD_, 'LineWidth',1,'DisplayName', 'Filtered sol.');
       
        ax1.XLim = [0 3];
        legend('box', 'off')
        ylabel('PSD');
        xlabel('Frequency');
        
        if ~isempty(folder)
            savefig(gcf,[folder,'PSD',suffix,'.fig']);
        end
    end
    if contains(fig_type,'m_abs_error')%-----------------------------------------------
        error       =   abs((V.p_xf - V.p_xf_true(P.t1KF:end)))/max(V.p_xf_true);
        t_error     =   S.t;
        box on;    hold on
         plot(t_error,error,...
             'DisplayName',['$\Delta t_\mathrm{analysis}=$ ',num2str(P.k_meas * P.dt)]) 

%              'DisplayName',['$\sigma_\mathrm{mic}=$ ',num2str(P.mic_sig_frac)])
%              'DisplayName',['$\sigma_\mathrm{frac}=$ ',num2str(P.sig)]) 
%              'DisplayName',['$N_\mathrm{mic}=$ ',num2str(P.N_mic)]) 
%              'DisplayName',['$m=$ ',num2str(P.m)])  
         legend('Orientation','horizontal','Location',...
            'northoutside','EdgeColor',[1 1 1]);
        xlabel('$t$')
        ylabel('$\frac{|p''_{\mathrm{f,ref}}-p''_{\mathrm{f,filter}}|}{\mathrm{max}(p''_{\mathrm{f,ref}})}$')

    end
    if contains(fig_type,'m_Cpp_error') %-----------------------------------------------       
        vec = [1,P.k_meas:P.k_meas:length(S.Cpp_tr)];
        t_Cpp_err   =   S.t(vec);
        Cpp_err     =   S.Cpp_tr(vec);        
        OUT.t       =   t_Cpp_err;
        OUT.Cpp_err	=   Cpp_err;          
    end
    if contains(fig_type, 'butterfly')   
        Aa      =   S.Aa_true;
        t_true 	=   S.t_true;
        % Compute the pressure response
        p_xf  	=   - P.sin_jpixf * Aa(P.N_m+1:2*P.N_m,:);
        % Perturb the solution at time t_p        
        t_      =   980 / P.dt; % Select some time here
        t_p     =   t_true;
        IC_p   	=   [Aa(1:P.N,t_) + 1e-6 * norm(Aa(1:P.N,t_),1); ...
                     Aa(P.N+1:end,t_)];
        [~,Aa_p]=   ode45(@(t,y)gov_eqns_adv_PE(t,y,P),t_p(t_:end),IC_p);

        Aa_b            =   Aa;
        Aa_b(:,t_:end)	=   Aa_p';
        p_xf_p          =   - P.sin_jpixf * Aa_b(P.N_m+1:2*P.N_m,:);

        
        % PLOT TIME EVOLUTION
        nexttile; plot(t_true, p_xf, t_p, p_xf_p)
        legend('Reference signal','Perturbed signal','Location',...
             'northoutside', 'Orientation', 'horizontal','Box', 'off')
        xlim([t_*P.dt, t_*P.dt+20])
        xlabel('$t$');     ylabel('$p''_\mathrm{f}$');

        % PLOT LOG DIFFERENCE
        diff = vecnorm(Aa - Aa_b,1,1);
        nexttile; plot(t_true,log(diff))
        xlim([t_* P.dt, t_ * P.dt + 50])
        hold on;
        axis square manual
        
        tA  =   t_ * P.dt+0.5; %x(1); % These values were selected manually first
        tB  =   tA + 10; %x(2); %


        xlim([tA, tA+80]); ylim([-10,7.5])
        xlabel('$t$');     ylabel('ln${||{\Delta\psi}||}$');

        poly    =   polyfit(t_true(tA/P.dt:tB/P.dt), ...
                            log(diff(tA/P.dt:tB/P.dt)), 1);
        f1      =   polyval(poly,t_true);
        plot(t_true,f1,'r--', 'LineWidth',2)
    end
end

%%





