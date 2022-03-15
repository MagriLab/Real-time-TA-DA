function [x_new,t_new] = fn_4thRunge_Kutta(t,x,dt,fn)
    % 4th order Runge kutta scheme:
    a = fn(t       , x         );
    b = fn(t + dt/2, x + a*dt/2);
    c = fn(t + dt/2, x + b*dt/2);
    d = fn(t + dt  , x + c*dt  );
    x_new = x + dt/6 * (a + 2*b + 2*c + d);
    t_new = t + dt;
end