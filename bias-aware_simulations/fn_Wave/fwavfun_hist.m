function [v] = fwavfun_hist(t,Vs,ts,j,J)
    try
        % Interpolates to find the wave function values at previous times:
        if t <= 0 || sum(ts) == 0
            % Assume 0 value before t = 0.
            v = 0;
        else
            % Interpolate to find previous values:
            st = max(1,j-J);
            ts = ts(st:j);
            Vs = Vs(st:j);
            if ts(end) == 0
                ts(end) = [];
                Vs(end) = [];
            end
            v = interp1(ts,Vs,t);
        end
    catch
        v = 0; 
        disp('no history provided to fwavfun_hist')
    end
end