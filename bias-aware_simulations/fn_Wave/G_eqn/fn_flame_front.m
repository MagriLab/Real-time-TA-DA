function [dxi_dt]  = fn_flame_front(varargin)%t,xi,ts,Xis,Gs,Hs,j,J,Mean,Geom)
    % Computes the Flame Front as in Eq. (2.8)
    if nargin == 10
        [t,xi,ts,Xis,Gs,Hs,j,J,Mean,Geom]=deal(varargin{:});
    elseif nargin == 5
        [t,xi,Sim,j,J]=deal(varargin{:});        
        Gs = Sim.Gs;
        Hs = Sim.Hs;
        ts = Sim.ts;
        Xis = Sim.Xis;
        Geom = Sim.Geom;
        Mean = Sim.Mean; 
    else
        error('wron number of inputs to fn_flame_front.m')
    end
    % Retreive
    Su = Mean.Su;
    
    % Wave functions
    [gt,~,gt_Tu] = fun_wavefun(t,ts,Xis,Gs,Hs,j,J,Mean,Geom);
    
    % Velocity definitions
    Ut      =   fn_flame_holder_vel(t,gt,gt_Tu,Mean,Geom);
    
    dxi_dr  =   fn_dxi_dr(xi,Ut,Mean,Geom);
    
    % Derivative of flame front wrt t:
    dxi_dt = Ut*ones(size(xi)) - Su*sqrt(1 + (dxi_dr).^2);
end