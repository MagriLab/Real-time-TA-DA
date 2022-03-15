function [U, ra] = advance_ESN(loop, varargin)
    ESN_P       =   varargin{1};
    if contains(loop, 'open')
        U     	=   varargin{2}; % washout
        [U,ra]	=   open_loop(U, ESN_P);
%         U       =   U(end,:);
    else
        ra   	=   varargin{2};
        U0   	=   varargin{3};
        N    	=   varargin{4};    % Number of integration timesteps
        [U, ra]	=   closed_loop(N, U0, ra, ESN_P);        
    end
    
end

% ======================================================================= %

function [U,ra] = open_loop(U0, ESN_P)
% ----------------------------------------------------------------------- %
%     Advances ESN in open-loop.
%         Args:
%             U: input bias time series (from data)
%         Returns:
%             ra: time series of augmented reservoir states
% ----------------------------------------------------------------------- %
%     U0(end,:)
%     U0          =   U0 + normrnd(0,ESN_P.std_noise);   
%     U0(end,:)
    N_units     =   ESN_P.N_units;
    N           =   size(U0,1);
    bias_out    =   1;
%     ra          =   zeros([1,N_units-1]);
    Wout    =   squeeze(ESN_P.Wout);
    U  =   U0;
    ra          =   horzcat(zeros([1,N_units-1]), bias_out);
    for i = 2:N
        ra = step(ra, U0(i-1,:), ESN_P);
        U(i,:)  =   ra * Wout;
    end

end

function [U, ra] = closed_loop(N, U0, ra, ESN_P)
% ----------------------------------------------------------------------- %
%     Advances ESN in closed-loop.
%         Args:
%             N: number of time steps
%             ra: initial reservoir state
%             Wout: output matrix
%         Returns:
%             U: time series of predicted bias
%             ra: final augmented reservoir state
% ----------------------------------------------------------------------- %
    U(:,1)  =   U0;
    Wout    =   squeeze(ESN_P.Wout);
    for i = 2:N
        ra      =   step(ra, U(:,i-1), ESN_P);
        U(:,i)  =   ra * Wout;
    end
end


function r_augmented = step(r_pre, U, ESN_P)
% ----------------------------------------------------------------------- %
%     Advances one ESN time step.
%         Args:
%             r_pre: reservoir state
%             U: input
%         Returns:
%             r_augmented: new augmented state (with bias_out appended)
% ----------------------------------------------------------------------- %
    % Extract hyperparameters and W matrices
    rho         =   ESN_P.hyperparameters(1);
    sigma_in    =   ESN_P.hyperparameters(2);
    bias_in     =   ESN_P.hyperparameters(end);
    bias_out    =   1;
    W           =   squeeze(ESN_P.W);
    Win         =   squeeze(ESN_P.Win);
    norm        =   ESN_P.norm;
    % Force U to be a row vector
    if size(U,1) ~= 1; U = U'; end
    % input is normalized and input bias added
    u_augmented =   horzcat(U./norm, bias_in); 
    % hyperparameters are explicit here
    r_pre       =   r_pre(1:end-1);
    r_post      =   tanh(sigma_in * u_augmented * Win + rho * r_pre * W);
    % output bias added
    r_augmented =   horzcat(r_post, bias_out);                              
end

