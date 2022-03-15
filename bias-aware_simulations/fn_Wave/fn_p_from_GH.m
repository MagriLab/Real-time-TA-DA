
function [ups,pps,rhops] = fn_p_from_GH(varargin)%t,x,ts,Gs,Hs,Mean,Geom,Acoustic_only)
    % Computes the modeshapes at t
    if nargin == 2        
        mi      =   varargin{1};    
        Sim     =   varargin{2};          
        Gs      =   Sim.Ensemble.Gs(mi,:);
        Hs      =   Sim.Ensemble.Hs(mi,:);
        ts      =   Sim.ts;
        index   =   find(Gs,1,'last');
        x_mic   =   Sim.Measurement.Mic_Pos;
        t       =   Sim.ts(index);
        Mean    =   Sim.Mean;
    elseif nargin == 7
        [t, x_mic, Mean, Geom, ts, Gs, Hs]     =   deal(varargin{:});
    else
        error('wrong number of inputs to fn_p_from_GH.m')
    end
    

    i  	=   0;
    [pps, ups, rhops]	=   deal(zeros(length(x_mic),1));
    for x = x_mic; i = i + 1;
        if x <= 0
            % Retrieve variables:
            c1      =   Mean.c1;
            u1      =   Mean.u1;
            rho1    =   Mean.rho1;
            Tu      =   Mean.Tu;
            try
                R_in    =   Mean.R_in;
            catch
                R_in    =   (1-Mean.M1)/(1+Mean.M1);
            end
            % Definition of f and wave functions
            g       =   @(t) fv_interp(t,ts,Gs);
            f       =   @(t) R_in * g(t - Tu);
            F       =   f(t - x/(c1 + u1));
            G       =   g(t + x/(c1 - u1));
            % Modeshapes
            pp      =   F + G;
            up      =   1/(rho1*c1)*(F - G);
            rhop    =   NaN;
        else
            % Retrieve variables:
            c2      =   Mean.c2;
            u2      =   Mean.u2;
            rho2    =   Mean.rho2;
            Td      =   Mean.Td;
            try
                R_out   =   Mean.R_out;
            catch
                R_out = -1;
            end
            % Definition of j and wave functions
            h       =   @(t) fv_interp(t,ts,Hs);
            j       =   @(t) R_out * h(t - Td);
            H       =   h(t - x/(c2 + u2));
            J       =   j(t + x/(c2 - u2));
            % Modeshapes
            pp      =   H + J;
            up      =   1/(rho2*c2)*(H - J);
            rhop    =   NaN; % Missing entropy waves
        end
        pps(i)      =   pp;
        ups(i)      =   up;
        rhops(i)    =   rhop;
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