function [Aa_m, Psi_a, py] = EnSRKF(Af_, py, Cee, Filter)
% Algorithm for the Ensemble Square-Root Kalman filter. 
% ----------------------------------------------------------------------- %
% Inputs:
%   -   Af_: forecast state matrix
%   -   py: observations 
%   -   Cee: observations covariance matrix
%   -   P: struct of parameter of the problem
% Outputs:
%   -   Aa_m: mean of the analysis ensemble matrix 
%   -   Psi_acomb: deviations of the analysis ensemble matrix
%   -   py: observations with the std applied
% ======================================================================= %
    % Retrieve variables
    M       =   Filter.M;
    m       =   Filter.m;
    % Ensemble mean
    psi_f_m	=   mean(Af_, 2); 
    % Create matrices of observations and mean
    Y  	=   py';
    Af_m	=   psi_f_m;         
    % Compute the analysis state mean value
    Psi_f   =   (Af_ - Af_m);
%         Psi_f   =   P.inflation * (Af_ - Af_m);
    if ~isreal(Psi_f)
        Psi_f   =   Af_ - Af_m + 1E-5;
        disp('Not real')
    end
    % Transform into observation space
    S       =   M * Psi_f; 
    W       =   (m - 1) * Cee + S * S.'; 
    % Check condition number
    try cond_W  =   cond(W);
        if cond_W > 1e15
            K       =   S.' * pinv(W);
        else
            K       =   S.'/W;
        end
    catch
        disp('W ILL CONDITIONED!')
        Psi_a   =   NaN(size(Psi_f));
        Aa_m    =   NaN(size(Af_m));
        return
    end
    % Compute the analysis state mean value
    Aa_m    =   Af_m +  Psi_f * K * (Y - M * Af_m);
    % Obtain analysis ensemble perturbations
    VEV    	=   K * S;
    [V,E]  	=   eig(VEV,'matrix');   V  =   real(V); E   =   real(E) ;
    IE      =   (eye(m)  - E);
    if min(IE,[],'all') > -1E-6
        IE  =   abs(IE);
    end
    Psi_a	=	Psi_f * V * (IE^(0.5)) * V.';
%     if ~isreal(Psi_a)
%         disp([num2str(min(eye(P.m)  - E))), ' changed to positive') 
%         Psi_a	=	Psi_f * V * (abs(eye(P.m)  - E)^(0.5)) * V.';
%     end        
% ----------------------------------------------------------------------- %
end