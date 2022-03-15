function [t_mic, p_mic, u_mic] = psi_to_pu(t_vec,psi_vec,P,varargin)
% ----------------------------------------------------------------------- %
    c_0         =   P.Mean.c_0;    
    sin_omj_mic =   sin(P.omega_j./c_0 .* P.x_mic)';
    cos_omj_mic	=   cos(P.omega_j./c_0 .* P.x_mic)';
    
    u_mic	=   cos_omj_mic * psi_vec(:,1:P.N_m)';
    p_mic	=   - sin_omj_mic * psi_vec(:,P.N_m+1:2*P.N_m)'; 
    
    if nargin == 3
        Fs_m    =   P.Fs_mic;
    else
        Fs_m    =   varargin{1};
    end
    t_mic 	=   t_vec(1):1/Fs_m:t_vec(end);
    
    % Interpolate u and p from the simulation to measurement time

    u_mic =  interp1(t_vec,u_mic',t_mic);
    p_mic =  interp1(t_vec,p_mic',t_mic);
end