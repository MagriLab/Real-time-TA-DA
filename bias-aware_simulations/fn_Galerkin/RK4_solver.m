function [t_next, psi_next] = RK4_solver(t,psi,dt,fun)
% Fourth order Runge–Kutta method for a general state vector and funciton.
% ------------------------------------------------------------------------ 
% Inputs:
%   - t: time of integration
%   - psi: psi(t) state vector at time t
%   - dt: time step
%   - fun: function of the governing equations
% Outputs:
%   - t_new: t+dt next time location
%   - psi_new: RK4 approximation of psi(t+dt)
% ------------------------------------------------------------------------
    % Fourth order RK coefficients
    K1  =   fun(t, psi); % Euler method
    K2  =   fun(t + dt/2, psi + K1*dt/2);
    K3  =   fun(t + dt/2, psi + K2*dt/2);
    K4  =   fun(t + dt, psi + K3*dt  );
    % Update the sate vector and time
    psi_next    =   psi + dt/6 * (K1 + 2*K2 + 2*K3 + K4);
    t_next      =   t + dt;
end