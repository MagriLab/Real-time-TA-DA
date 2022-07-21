
function portraits(p_cell, T_cell, row)

colors     	=   {' #1a5276 ',' #a3e4d7 ',' #f39c12 '};
markersizes	=   [24,16,12];
lines 	=   ['1','0.5','0.5'];
names = {'Truth','Filtered sol.','Unfiltered sol.'};
% figure('Units','normalized','OuterPosition',[0 0.06 0.8 0.8]); 
% tiledlayout(1,2,'Padding','compact','TileSpacing','compact')
for j = 1:length(p_cell)
    p_f = p_cell{j};
    T = T_cell{j};
    c = colors{j};
    ms = markersizes(j);
    l = lines(j);
    %% START ANALYSIS
    Fs  =   1/(T(2)-T(1));           % Sampling frequency
    N   =   length(T);
    P   =   p_f - mean(p_f);
    TP  =   T;
    %  Get rid of transient
    Nmax = 50000;
    if (N > Nmax) %&& (strcmp(truncate,'on'))
%         disp(['Reducing datapoints from ',...
%             num2str(N),' to ',num2str(Nmax),' to discard transient'])
        N   =   Nmax; %x=x(end-N+1:end); t=t(end-N+1:end);
        n   =   length(P)- N + 1;
        X   =   zeros(round(N),1);
        T   =   zeros(round(N),1);
        % Increase interval by neglecting every 10th point
        i   =   1;

        while(n <=  length(P))
            X(i)    =   P(n);
            T(i)    =   TP(n);
            n       =   n + 1;
            i       =   i + 1;
        end
    Fs  =   1/(T(2)-T(1));           % Sampling frequency
    % size(X)
    else
       X = P;
       T = TP;
    end

    if j==1
    %% Plot average mutual information and get optimum time delay
    lag_max     =   0.5*Fs; %Maximum lag
    % Average mutual information and lag positions from 'ARFM_ami.m' function
    [v,L]   =   ARFM_ami(X,lag_max);
    % The first local minimum is at the optimal time delay
    [~,tau_AMI]   =   findpeaks(-v);
    Delay               =   tau_AMI(1) + 0.2;
    zeta                 =   floor(Fs*Delay/1000);    % Transform to seconds 


    %% Obtaining embedding dimension from a False Nearest Neighbor method
    % Set typical input values for algorithm
    dmax=10; Rtol=50; Atol=2;
    % Extract the percentage value of false nearest neighbors from ARFM_fnn.m
    X_ = P(end-2000:end);
    T_ = TP(end-2000:end);
    [FNN]=ARFM_fnn(X_,T_,Fs,dmax,zeta,Rtol,Atol);

    fnn_zero = find(FNN==0,1,'first');
    d = fnn_zero + 1;

    end
    %% Phase space reconstruction
    % Find length of the delay vectors (M) using the values of number of data points in the signal (N),
    % minimum embedding dimension, (d) and the optimum time delay (tau) as M=N-(d-1)*tau
    M   =   length(X) - (d - 1) * zeta;
    % Find number of delayed vectors of a signal from the given values of embedding dimension and delay
    [Y]     =   ARFM_delay_vec(X,zeta,d,M);  

    ax = nexttile(row+1); hold on
    if size(Y,2) == 2
        plot(ax,Y(:,1),Y(:,2),'k','linewidth',l,'color',c); grid on;
        title('2D Phase Portrait')
    else
        plot3(ax,Y(:,1),Y(:,2),Y(:,3),'k','linewidth',1.5,'color',c); grid on;
        zlabel(ax,'$p''(t+2\zeta)$')
    end
    xlabel(ax,'$p''(t)$')
    ylabel(ax,'$p''(t+\zeta)$')
    axis(ax,'square');
    box(ax, 'on')
    view(ax, [180 0])


    %% Plot First Return Map
    [pks,locs] = findpeaks(X,'MinPeakDistance',10,'Minpeakheight',0);
    M=max(pks)+2*max(pks);


    ax = nexttile(row+2); hold on;
    if j==1
        plot(ax,0:1/Fs:10000,0:1/Fs:10000,'linewidth',0.5,'Color', 'k','HandleVisibility','off'); 
        grid on;  hold on;
    end
    plot(ax,pks(1:end-1),pks(2:end),'.','Color',c,'MarkerSize',ms,'DisplayName',names{j});
    % title('First return map')
    xlabel(ax,'$p''_\mathrm{max} (i)$')
    ylabel(ax,'$p''_\mathrm{max} (i+1)$')
    xlim(ax,[0,max(M,0.2)]);
    ylim(ax,[0,max(M,0.2)]);
    axis(ax,'square');
        legend('Location', 'eastoutside','NumColumns', 1,'EdgeColor','None')

    
end
end