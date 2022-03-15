function [ts,Gs,Hs,Xis,Qs] = fn_time_marching(t_end,dt,Geom,Mean)
    % Function to time march a simulation:
    
    % Retreive variables:
    Tau =   Mean.Tau;
    Tu  =   Mean.Tu;
    N   =   Geom.N;
    
    % Time marching paramters:
    Nt  =   ceil(t_end/dt) + 1;
    J   =   ceil((Tau + Tu)/dt) + 10;    
    
    % Initial conditions:
    xi  =   Mean.xibar;
    t   =   0;

    % Prelocate
    [Gs, Hs, Qs, ts] = deal(zeros(Nt,1));
    Xis = zeros(N,Nt);
    for j = 1:Nt
        if j == 1
            fprintf('\t\tProgress: 0%%,')
        elseif mod(j,round(Nt/4)) == 0
            if round(j/Nt,2)*100 < 100
                fprintf([' ',num2str(round(j/Nt,2)*100),'%%, '])
            end
        end
        idx = max(j-10*J,1);
        % Compute g(t) and h(t)
%         [gt,ht,~,Qt] = fun_wavefun(t,ts,Xis,Gs,Hs,j,J,Mean,Geom);
        [gt,ht,~,Qt] = fun_wavefun(t,...
                        ts(idx:j),Xis(:,idx:j),Gs(idx:j),Hs(idx:j),...
                        j-idx+1,J,Mean,Geom);
        Gs(j) = gt;
        Hs(j) = ht;
        Qs(j) = Qt;
        ts(j) = t;
        Xis(:,j) = xi;
        if Mean.Qbar == 0
            t = t + dt;
        elseif ~isfield(Mean, 'Heat_law') || strcmp(Mean.Heat_law,'G_eqn')
            % Compute new xi:
            fun = @(t,xi) fn_flame_front(t,xi,ts,Xis,Gs,Hs,j,J,Mean,Geom);
            % Time march
            [xi,t] = fn_4thRunge_Kutta(t,xi,dt,fun);
            if max(abs(xi)) > 1e4
                if Geom.Speaker
                    error(['Solution not found. ',...
                           'Change operating conditions or forcing.'])
                else
                    error(['Solution not found.',...
                           'Change operating conditions.'])
                end
            end
        else
            t = t + dt;
        end
        
    end
    fprintf('100%%.\n')
end