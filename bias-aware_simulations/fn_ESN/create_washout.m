
function [U_washout, t_washout] = create_washout(dt, N_washout, t_f, p_f, t_t, p_t)
    % Define washout time
    t_w_start   =  	t_t(end) - dt*double(N_washout);
    t_w_end     =	t_t(end);
    t_washout   = 	t_w_start:dt:t_w_end;
    % Evaluate pressure at washout timesteps
    p_t         =   interp1(t_t, p_t', t_washout, 'spline');
    p_f         =   interp1(t_f, p_f', t_washout, 'spline');
    % Bias
    U_washout   =   p_t - p_f;
    U_washout   =   U_washout(end-N_washout+1:end,:);
    t_washout   =   t_washout(end-N_washout+1:end);
end