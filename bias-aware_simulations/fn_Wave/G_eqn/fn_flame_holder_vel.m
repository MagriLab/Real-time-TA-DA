function [Ut] = fn_flame_holder_vel(t,gt,gt_Tu,Mean,Geom)
    % Computes the velocity at the flame holder
    A1      =   Geom.A1;    u1      =   Mean.u1;
    c1      =   Mean.c1;    M1      =   Mean.M1;
    rho1    =   Mean.rho1;  R_in    =   Mean.R_in;
    try if Geom.Speaker
            Tu    = Mean.Tu;
            V     = Mean.Force.Vu;
            omega = Mean.Force.omega_u;
            F = V*c1/(A1*(1 + M1))*cos(omega*(t - Tu/2));
        else; F = 0;
        end
    catch; F = 0; end
    ft = R_in * gt_Tu + F;
    % Velocity U(t)
    Ut = u1 + 1/(rho1*c1) * (ft - gt); 
end