%% Modeshapes
function [up,pp,rhop] = fn_modehshapes(varargin)%t,x,ts,Gs,Hs,Mean,Geom,Acoustic_only)
    % Computes the modeshapes at t, x
    if nargin == 8
        [t,x,ts,Gs,Hs,Mean,Geom,Acoustic_only] = deal(varargin{:});
    elseif nargin == 2
        Sim     =   varargin{1};
        Mean    =   Sim.Mean;        
        index   =   find(Sim.Gs(1,:),1,'last');
        Gs      =   Sim.Gs;
        Hs      =   Sim.Hs;
        ts      =   Sim.ts;
        t       =   ts(index);        
        x       =   varargin{2};
    else
        error('wrong number of inputs to fn_modehshapes.m')
    end
    % =================================================================== %
    if x <= 0
        % Retrieve variables:
        c1      =   Mean.c1;
        u1      =   Mean.u1;
        rho1    =   Mean.rho1;
        Tu      =   Mean.Tu;
        R_in    =   Mean.R_in;
        % Definition of f and wave functions
        g   =   @(t) fv_interp(t,ts,Gs);
        f   =   @(t) R_in * g(t - Tu);
        F   =   f(t - x/(c1 + u1));
        G   =   g(t + x/(c1 - u1));
        % Modeshapes
        pp      =   F + G;
        up      =   1/(rho1*c1)*(F - G);
%         rhop    =   rho1 + 1/(c1^2) * (F + G);
        rhop    =   NaN;
    else
        % Retrieve variables:
        c2      =   Mean.c2;
        u2      =   Mean.u2;
        rho2    =   Mean.rho2;
        Td      =   Mean.Td;
        R_out   =   Mean.R_out;
        % Definition of j and wave functions
        h   =   @(t) fv_interp(t,ts,Hs);
        j   =   @(t) R_out * h(t - Td);
        H   =   h(t - x/(c2 + u2));
        J   =   j(t + x/(c2 - u2));
        % Modeshapes
        pp  = H + J;
        up  = 1/(rho2*c2)*(H - J);
        rhop    = NaN; % Missing entropy waves
    end
end

%% Function definitions:
function v = fv_interp(t,ts,Vs)
    % Interpolate
    if t < 0
        v = 0;
    else
        v =  interp1(ts,Vs,t);
    end
end