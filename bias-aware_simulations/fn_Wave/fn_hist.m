function [v] = fn_hist(t,Vs,ts,J)
% Funtion that interpolates to find the values of a quantity Vs at a
% previous time t.
% ----------------------------------------------------------------------- %
    try
        if t <= 0
            % Assume 0 value before t = 0 or mean flame front.
            if min(size(Vs)) == 1
                v = 0;
            else
                v = Mean.xibar;
            end
            
        else
            % Interpolate to find previous values:
            st = max(1,length(ts)-J);
            ts = ts(st:end);
            Vs = Vs(:,st:end);
            if ts(end) == 0
                ts(end) = [];
                Vs(:,end) = [];
            end
            v = interp1(ts,Vs.',t).';
        end
    catch
        v = 0; 
        disp('no history provided to fwavfun_hist')
    end
end