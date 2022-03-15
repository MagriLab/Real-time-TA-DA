function omega_j = fn_obtain_omega(varargin)%N_m, L_0, x_f, c1_0, c2_0)
    
    % Note: the subscript _0 indicates mean flow property    
    if nargin == 1
        if isstruct(varargin{1})
            N_m     =   varargin{1}.N_m;
            try 
                c1_0    =   varargin{1}.c1_0;
                c2_0    =   varargin{1}.c2_0;
            catch
                c1_0    =   varargin{1}.Mean.c_0;
                c2_0    =   c1_0;
            end
                
            L_0     =   varargin{1}.L_0;
            x_f     =   varargin{1}.x_f;
        else
            N_m     =   varargin{1};
            c1_0    =   300;
            c2_0    =   350;
            L_0     =   1;
            x_f     =   0.25;
        end
    elseif nargin == 5
        [N_m, L_0, x_f, c1_0, c2_0] = deal(varargin{:});
    elseif nargin == 4
        [N_m, L_0, x_f, c1_0] = deal(varargin{:});
        c2_0    =   c1_0;
    end
    % ------------------------------------------------------------------- %
    % Relation between alpha_m and beta_m, jump condition
    fun     =   @(omega) c2_0 * sin(omega * x_f / c1_0) * ...
                        cos(omega * (L_0 - x_f) / c2_0) + ...
                        c1_0 * cos(omega * x_f / c1_0) * ...
                        sin(omega * (L_0 - x_f) / c2_0);
    
    % Weighted averaged mean density
    c_0 =   (1 - x_f/L_0) * c2_0 + x_f/L_0 * c1_0;
    
    
    % Find roots with approximate omega as initial guess
    omega_j =   zeros(1, N_m);
    IC      =   zeros(1, N_m);
    for j = 1:N_m
        IC(j)     	=   j * pi * c_0 / L_0;
        if fun(IC(j)) < 1e-6
            omega_j(j)      =   IC(j);
        else
            [omega_j(j),~]  =   fsolve(fun, IC(j));
        end
    end
%     figure; hold on;
%     plot(omega_j, omega_j*0, '+','MarkerSize', 10);
%     plot(IC, IC*0, 'x', 'MarkerSize', 10, 'color',[0.4940, 0.1840, 0.5560])
%     legend('$\omega_j$','Approx. $\omega_j$')
end