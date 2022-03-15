function [xi] = fn_flame_front_hist(t,Xis,ts,j,J,Mean)
    % Interpolates to find the flame front function values at previous
    % times:
    if t <= 0 || sum(ts) == 0
        % At times t< 0 ensure that the flame front is given by the mean.
        xi = Mean.xibar;
    else
        % Interpolate to find previous values:
        st = max(1,j-J);
        ts = ts(st:j);
        Xis = Xis(:,st:j);
        if ts(end) == 0 && j > 1
            ts(end) = [];
            Xis(:,end) = [];
        end
        xi = interp1(ts,Xis.',t).';
    end
end