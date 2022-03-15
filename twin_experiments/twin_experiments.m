function twin_experiments(choose_fig)
folder = dir('*/data_twin/'); folder = folder.folder;
save_  = true;
choose_fig = choose_fig - 3;
%% =================== CREATE AND SAVE or LOAD DATA ==================== %%
% plot_settings
timeseries_figs	=   [5, 9, 13, 16]; params_figs    	=   [7,18,19];
error_figs    	=   [6,10,15];      inflation_figs 	=   [11, 12];
butterfly_fig   =   14;             PSD_fig      	=   17;
file_name       =   [folder,'/Fig',num2str(choose_fig+3),'.mat'];
if ~isfolder(folder)
   mkdir(folder);   disp('Data folder created')
end
if isfile(file_name)
    disp('Loading data...')
    load(file_name); disp('done')
else % Create data
    disp('Creating data...')
    [ALL_S,ALL_V,ALL_P] = DA_paper_data(choose_fig); disp('done')
    if save_
        disp('Saving data...')
        save(file_name, 'ALL_S','ALL_V', 'ALL_P','-v7.3'); disp('done')
    end
end
%% ===================== PLOT THE SELECTED FIGURE ====================== %%   
disp('Plotting...')
% ------------------------------------------ Timeseries non-chaotic figures
if any(choose_fig == timeseries_figs)
    figure(11);
    tiledlayout(2,2,'TileSpacing','Compact');
    for i = 1:length(ALL_S)
        nexttile(i);   Figures(ALL_S{i},ALL_V{i},ALL_P{i},'timeseries',[]);
        YL= get(gca, 'YLim');
        set(gca,  'YLim', [YL(1), YL(2)*3.5]);
    end
    figure(12); hold on;
    tiledlayout(2,2,'TileSpacing','Compact');
    for i = 1:length(ALL_S)
        nexttile;	Figures(ALL_S{i},ALL_V{i},ALL_P{i},'zoom',[]);
    end
end   
% --------------------------------------------- Parameter evolution figures 
if any(choose_fig == params_figs)
    figure; tiledlayout(2,2,'TileSpacing','Compact');
    for i = 1:length(ALL_S)
        nexttile;   Figures(ALL_S{i},ALL_V{i},ALL_P{i},'params',[]);
    end
end
% ---------------------------------------- Increase, inflate,reject figures 
if any(choose_fig == inflation_figs)
    figure; tiledlayout(2,2);
    mi2 = 0;
    for vi = [1,3]
        for mi = 1:4
            mi2     =   mi2 + 1;
            if mi2 > size(ALL_S,2); break; end
            S = ALL_S{mi2};  P = ALL_P{mi2};
            nexttile(vi); hold on;
            title(['$\rho =$', num2str(P.inflation)]) 
            Figures(S,[],P,['params_2_', num2str(vi+1)],[]);
            title(['$\rho =$', num2str(P.inflation)]) 
        end
    end
end
% ------------------------------------------ Power spectral density figures 
if any(choose_fig == PSD_fig)
    figure(21); tiledlayout(2,2,'TileSpacing','Compact');
    tl_order = [1,3;2,4];
    for i = 1:length(ALL_S); j = 1;
        V = ALL_V{i}; P = ALL_P{i};
        init    =   [900, 1200];    w   =   init(2)-init(1);
        for t1 = init
            nexttile(tl_order(i,j));hold on            
            [f,PSD_] =   fun_PSD(P.dt,V.p_xf_true(t1/P.dt+1:(t1+w)/P.dt));
            plot(f, PSD_, 'LineWidth',4,'DisplayName', 'True sol.');
            [f,PSD_] =   fun_PSD(P.dt,V.p_xf_comp((t1-P.t_minKF)/P.dt+1:...
                                                  (t1-P.t_minKF+w)/P.dt));
            plot(f, PSD_,'LineWidth',0.8,'DisplayName', 'Unfiltered sol.');
            [f,PSD_] =   fun_PSD(P.dt,V.p_xf((t1-P.t_minKF)/P.dt+1:...
                                             (t1-P.t_minKF+w)/P.dt));
            plot(f, PSD_, 'LineWidth',1,'DisplayName', 'Filtered sol.');
            xlim([0,3]); ylim([0, 1.2])
            legend('box', 'off');  ylabel('PSD'); xlabel('Frequency');
            text(0.5,1.1,['$t\in$[',num2str(t1),',',num2str(t1+w),'$]$'])
            j=2;
        end
    end
end
% ----------------------------------------------- Lyapunox exponent figures 
if choose_fig == butterfly_fig
    figure; tiledlayout(2,2)
    Figures(ALL_S{1},[],ALL_P{1},'butterfly',[])    
end
% --------------------------------------------------- Error metrics figures
if any(choose_fig == error_figs)
    figure; tiledlayout(2,2,'TileSpacing','Compact');
    ax	=   nexttile(1); hold(ax,'on');
    legend('Orientation','horizontal','Location','northoutside',...
        'EdgeColor',[1 1 1]);
    xlabel('$t$'); 
    ylabel(['$\frac{|p''_{\mathrm{f,ref}}-p''',...
            '_{\mathrm{f,filter}}|}{\mathrm{max}(p''_{\mathrm{f,ref}})}$'])
    ALL_Cpp_err = []; ALL_Cpp_t = [];
    %%%%%% LEFT HAND SIDE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 1:3; S = ALL_S{i};    V = ALL_V{i};      P = ALL_P{i};     
        labl =   '$\sigma_\mathrm{frac}= $';
        if choose_fig == 6
            labl    =   [labl, num2str(P.sig)];
        elseif choose_fig == 10
            labl    =   [labl, num2str(P.mic_sig_frac)];
        elseif choose_fig == 15
            labl    =   [labl, num2str(P.k_meas * P.dt)]; 
        end
        % Store covariance trace
        ALL_Cpp_t	=   [ALL_Cpp_t, S.t(:)];
        ALL_Cpp_err	=   [ALL_Cpp_err, S.Cpp_tr(:)];
        % Compute relative error
        error	=   abs((V.p_xf-V.p_xf_true(P.t1KF:end)))/max(V.p_xf_true);
        t_err	=   S.t;
        % PLOT RELATIVE ERROR
        ax	=   nexttile(1); hold(ax,'on');	box on; hold on
        plot(t_err,error,'DisplayName',labl)
    end
    % PLOT COVARIANCE TRACE       
    ax	=   nexttile(3); hold(ax,'on')
    nz  =   find(~(ALL_Cpp_err(:,1)==0));
    tt  =   ALL_Cpp_t(nz,1);
    yy  =   ALL_Cpp_err(nz,:);
    bar(tt, yy, 'grouped','edgecolor','w')
    set(gca,'YMinorTick','on','YScale','log'); xlim([900,P.t_maxKF])     
    xlabel('$t$'); 
    ylabel('tr($\mathrm{\mathbf{C_{\psi\psi}^\mathrm{f}}}$)')
    %%%%%% RIGHT HAND SIDE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    if choose_fig ~= 15
        ax	=   nexttile(2); hold(ax,'on')
        legend('Orientation','horizontal','Location','northoutside',...
            'EdgeColor',[1 1 1]);
        xlabel('$t$'); ylabel(['$\frac{|p''_{\mathrm{f,ref}}-p',...
            '_{\mathrm{f,filter}}|}{\mathrm{max}(p''_{\mathrm{f,ref}})}$'])
        ALL_Cpp_err = [];   ALL_Cpp_t = [];
        for i = 4:6; S = ALL_S{i};    V = ALL_V{i};      P = ALL_P{i};       
            if choose_fig == 6
                labl    =   ['$m=$ ',num2str(P.m)];
            elseif choose_fig == 10
                labl    =   ['$\Delta t_\mathrm{analysis}=$ ',...
                            num2str(P.k_meas * P.dt)]; 
            end
            % Store covariance trace
            ALL_Cpp_t	=   [ALL_Cpp_t, S.t(:)];
            ALL_Cpp_err	=   [ALL_Cpp_err, S.Cpp_tr(:)];
            % Compute relative error
            error	=   abs((V.p_xf - V.p_xf_true(P.t1KF:end)))/...
                            max(V.p_xf_true);
            t_err   =   S.t;
            % PLOT RELATIVE ERROR
            ax	=   nexttile(2); hold(ax,'on');	box on; hold on
            plot(t_err,error,'DisplayName',labl)
        end     
        % PLOT COVARIANCE TRACE          
        ax	=   nexttile(4); hold(ax,'on')
        nz  =   find(~(ALL_Cpp_err(:,1)==0));
        tt  =   ALL_Cpp_t(nz,1);
        yy  =   ALL_Cpp_err(nz,:);
        bar(tt, yy, 'grouped')        
        set(gca,'YMinorTick','on','YScale','log'); xlim([900,P.t_maxKF])     
        xlabel('$t$'); 
        ylabel('tr($\mathrm{\mathbf{C_{\psi\psi}^\mathrm{f}}}$)')
    end
end
fprintf('done\n')
