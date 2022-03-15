function fn_plot_mean_std(t, y, name)

    y = squeeze(y);
        
    mean_   =   mean(y)';
    std_    =   std(y)';
    std_p   =   mean_+ std_;
    std_m   =   mean_- std_;
    l   =   plot(t(:),mean_,'LineWidth',1.4,'DisplayName', name); 
    c   =   get(l,'Color');
    patch([t(:); flipud(t(:))],[std_m; flipud(std_p)],c,...
        'FaceAlpha',0.2, 'EdgeColor','none','DisplayName',['std(',name,')']); 
    
end