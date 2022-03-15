function fn_plot_results(Filter, Sim_t, Truth, t_KF, Aa_KF, Obs, t_bias, bias, t_o, U_o)
    folder = dir(fullfile('./','*/higher_fidelity_simulations.m')); 
    folder = [folder(1).folder,'\Figures\'];
    if ~isfolder(folder); mkdir(folder); end
    prompt = "\n Do you want to save the figures? Y/N [Y]: ";
    saving = upper(input(prompt,"s"));
    if isempty(saving) || saving == 'Y'
        saving = true;
        prename = 'Matlab_';
        if Filter.E_State; prename(end+1) = 'S'; end
        if Filter.E_Params; prename(end+1) = 'P'; end
        prename = [prename, 'E_bias_'];
    else;  saving = false;
    end
    %% ================== TIMESERIES, PARAMETERS & PSD ================== %
    Aa 	=   squeeze(mean(Aa_KF, 1));
    p	=   - Filter.sin_omj_mic * Aa(Filter.N_m+1:2*Filter.N_m,:); 
    colors     	=   {' #1a5276 ',' #a3e4d7 ',' #f39c12 '};
    mic     =   1; N_m = Filter.N_m;
    xlims   =   [Filter.t_start-0.04, Filter.t_stop+0.1];
    figure('Units','normalized','OuterPosition',[0 0.06 1 0.95]); 
    tiledlayout(2,2,'Padding','compact','TileSpacing','compact')
    f1=nexttile; hold on%--------------------------------------- TIMESERIES
    l=plot(Sim_t.t_mic,Sim_t.p_mic(mic,:),'LineWidth',3,'color','k'); 
    l.Color(4)=0.15;
    l=plot(t_KF, p(mic,:),'--','color','k'); l.Color(4)=0.4;
    if ~isempty(bias)
        bias_interp	=   interp1(t_bias, bias', t_KF)';
        bias_interp(isnan(bias_interp))=0;
        p_unbias 	=   p + bias_interp;
        plot(t_KF, p_unbias(mic,:),'-','color','k','LineWidth', 0.5);
    end
    xlabel('$\tilde{t}$ [s]'); ylabel('$p''_\mathrm{mic_1}$ [Pa]')
    xlim(xlims); 
    ylim([min(Sim_t.p_mic(mic,:))*1.1, max(Sim_t.p_mic(mic,:)) * 5])
    plot(Truth.t_mic, Obs(mic,:),'o','markersize',4,...
        'MarkerFaceColor','r','MarkerEdgeColor','k')
    if ~isempty(bias)
        l=legend('Truth', 'Biased filter sol.', ...
            'Unbiased filter sol.','Analysis step');
    else; l=legend('Truth', 'Filter sol.', 'Analysis step');
    end
    l.Location = 'NorthOutside'; l.Orientation = 'horizontal'; 
    l.EdgeColor = 'None';    lims = get(gca, 'YLim');
    if Filter.E_Params;  plot(Truth.t_mic(Filter.init_SPE_i)*[1,1], ...
                             lims, 'k--','HandleVisibility','off');
    end
    plot(Truth.t_mic(Filter.init_BE_i)*[1,1], lims, ...
        'k--','HandleVisibility','off');
    plot(Filter.t_stop*[1,1], lims, 'k--','HandleVisibility','off');
    f2=nexttile; hold on%--------------------------------------------- ZOOM
    xlim([1.4 1.44]); ax=gca; ax.FontSize = 28;
    co = copyobj(f1.Children,f2); 
    co(1).MarkerSize = 10; co(2).LineWidth = 1.5;
    co(3).LineWidth = 2;co(4).LineWidth = 10;
    nexttile; hold on; xlabel('$t$ [s]')%--------------------------- PARAMS
    if Filter.E_Params == 1
        fn_plot_mean_std(t_KF, Aa_KF(:,end,:), '$\beta (t)$')
        ylabel(['$\tilde{\beta}$ $\left[{\mathrm{W}}',...
            '{\frac{\mathrm{s}^{1/2}}{\mathrm{m}^{5/2}}}\right]$'])
    elseif Filter.E_Params == 2
        fn_plot_mean_std(t_KF, Aa_KF(:,end,:)/Filter.tau,...
                        '$\tau (t)/\tau^\mathrm{true}$')
        fn_plot_mean_std(t_KF, Aa_KF(:,end-1,:)/Filter.beta, ...
                        '$\beta (t)/\beta^\mathrm{true}$')
        l=legend(); l.Location = 'EastOutside'; 
        l.Orientation = 'horizontal'; l.EdgeColor = 'None'; 
    end
    lims = get(gca, 'YLim'); t = Truth.t_mic;
    plot(t(Filter.init_SPE_i)*[1,1],lims, 'k--','HandleVisibility','off');
    plot(t(Filter.init_BE_i)*[1,1],lims, 'k--','HandleVisibility','off');
    plot(Filter.t_stop*[1,1], lims, 'k--','HandleVisibility','off');
    ylim(lims); xlim(xlims)
    try  true_p1	=	- Filter.sin_omj_mic * Sim_t.psi(:,N_m+1:2*N_m)';
    catch; true_p1	=   Sim_t.p_mic;
    end
    true_p1     =   interp1(Sim_t.t_mic, true_p1', t_KF)';
    actual_bias =   true_p1(mic,:) -  p(mic,:);
    %
    f4=nexttile(4); hold on %------------------------------------------ PSD
    %
    t1          =   1.2; [~,loc1]   =   min(abs(t_KF - t1));
    t2          =   1.6; [~,loc2]   =   min(abs(t_KF - t2));
    dt          =   1/Sim_t.Measurement.Fs; 
    [f,PSD_]    =   fun_PSD(dt,true_p1(mic,loc1:loc2)');
    plot(f, PSD_,'color', colors{1},'LineWidth',3); 
    [f,PSD_]    =   fun_PSD(dt,p(mic,loc1:loc2)');
    plot(f, PSD_,'color', colors{2},'LineWidth',1.5);
    [f,PSD_]    =   fun_PSD(dt,p_unbias(mic,loc1:loc2)');
    plot(f, PSD_,'color', colors{3},'LineWidth',1.5);
%     f_up    =    343/(1.18+0.74)/2;     f_down  =    589/(1.18+0.74)/2;
    ylabel('PSD'); xlabel('Frequency [Hz]'); 
    xlim([0,300]); yl = get(gca,'YLim');
    text(5,yl(2)*0.9, ['$\tilde{t} \in [',num2str(t1,2),', ',...
         num2str(t2,2),']$ s'],'HorizontalAlignment','left')
    l=legend('Truth','Biased filter sol.','Unbiased filter sol.');
    l.EdgeColor = 'none';

    if saving
        fn_save_pdf_fig(folder, [prename, 'timeseries'])
    end
    %% ============= PHASE PORTRAIT AND FIRST RETURN MAP ================ %
    p_cell = {true_p1(mic,loc1:loc2),p(mic,loc1:loc2),... 
              p_unbias(mic,loc1:loc2)};
    t_cell = {t_KF(loc1:loc2), t_KF(loc1:loc2), t_KF(loc1:loc2)};
    plot_portraits(p_cell,t_cell)
    if saving
        fn_save_pdf_fig(folder, [prename, 'portrait'])
    end
    %% ==================== ACTUAL BIAS vs ESN BIAS ===================== %
    figure('Units','normalized','OuterPosition',[0 0.06 1 0.45]); 
    tiledlayout(1,2,'Padding','compact','TileSpacing','compact')
    f1=nexttile;hold on %--------------------------------------------------
    l = plot(t_KF, actual_bias, 'LineWidth', 2, 'color', '#4DBEEE'); 
    l.Color(4) = 0.8;
    l2 = plot(t_bias, squeeze(bias(mic,:)),'--', 'color',...
               '#7E2F8E','LineWidth', 1.5);
    plot(t_o, U_o(:,1), '-.', 'Color', l2.Color, 'LineWidth', 1.5)
    xlim([1 1.198]); ylim([min(actual_bias) max(actual_bias)]*1.02)
    xlabel('$\tilde{t}$ [s]')
    ylabel('BIAS $= \Delta \tilde{p}_\mathrm{mic_1}''$ [Pa]')

    reinit = interp1(t_KF,actual_bias ,Truth.t_mic(Filter.init_BE_i:end))';
    plot(Truth.t_mic(Filter.init_BE_i:end), reinit, 'o',...
        'markersize',6,'MarkerEdgeColor', l2.Color)
    f2=nexttile;hold on %--------------------------------------------------
    xlim([Truth.t_mic(end), Truth.t_mic(end)+0.20])
    copyobj(f1.Children, f2)
    xlabel('$\tilde{t}$ [s]'); f2.YTickLabel = '';
    legend('Actual bias', 'Closed-loop ESN','Open-loop ESN',...
           'Analysis step','Location', 'northoutside', 'edgecolor', ...
           'none', 'orientation','horizontal')
    if saving
        fn_save_pdf_fig(folder, [prename, 'ESN'])
    end
end
