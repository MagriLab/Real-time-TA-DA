function [bias, t_bias] = fn_plot_bias(Sim_t, Sim_f, varargin)
    
    t_HOM = Sim_t.t_mic;
    p_HOM = Sim_t.p_mic;
    t_LOM = Sim_f.t_mic;
    p_LOM = Sim_f.p_mic;
    
    if size(p_LOM,1)<size(p_LOM,2)
        p_LOM = p_LOM';
    end
    if size(p_HOM,1)<size(p_HOM,2)
        p_HOM = p_HOM';
    end
        
    t1      =   1/(t_HOM(2)-t_HOM(1));
    if length(t_HOM) < 2*t1
        t1 = 1;
    end
    try
        bias    =   (p_HOM(t1:end,:) - p_LOM(t1:end,:))';
    catch
        p_HOM   =   interp1(t_HOM, p_HOM, t_LOM);
        t_HOM   =   t_LOM;
        bias    =   (p_HOM(t1:end,:) - p_LOM(t1:end,:))';
    end 
    t_bias 	=   t_HOM(t1:end);

    if nargin == 3
        filename = varargin{1};
        if ~isfolder('ESN_bias_data')
            mkdir('ESN_bias_data')
        end
        if ~isfile([filename,'.mat']); save([filename,'.mat'],'bias')
        %else; fprintf('\nData not saved. filename already exists\n')
        end
    end
    
    tzoom   =   0.05/(t_HOM(2)-t_HOM(1));       
    figure; tiledlayout(2,2)
    c = [0 0.4470 0.7410;0.4660 0.6740 0.1880;0.6350 0.0780 0.1840];
    nexttile; hold on
    plot(t_HOM(t1:end),p_HOM(t1:end,1), 'color', c(1,:))
    ylabel('$p_\mathrm{mic_1}$ [Pa]'); xlim([t_bias(1), t_bias(end)])
    plot(t_LOM(t1:end),p_LOM(t1:end,1), 'color', c(2,:))
    legend('truth', 'LOM')
    nexttile; hold on 
    plot(t_HOM(end-tzoom:end),p_HOM(end-tzoom:end,1), 'color', c(1,:))      
    plot(t_LOM(end-tzoom:end),p_LOM(end-tzoom:end,1), 'color', c(2,:))
    xlim([t_LOM(end-tzoom) t_LOM(end)])
    legend('truth', 'LOM')
    nexttile; plot(t_bias, bias(1,:), 'color', c(3,:))    
    xlabel('time [s]')
    ylabel('BIAS $= \Delta p_\mathrm{mic_1}$ [Pa]'); xlim([t_bias(1), t_bias(end)])
    nexttile; plot(t_bias(end-tzoom:end), bias(1,end-tzoom:end), 'color', c(3,:))
    xlabel('time [s]')
end