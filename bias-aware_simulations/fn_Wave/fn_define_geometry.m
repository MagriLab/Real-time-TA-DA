function [Geom] = fn_define_geometry(order)
    % Provides the geometric definitions of the combustor:
    
    Lu = 1.18;     % Length of upstream working section (m)
    Lb = 0.74;     % Length of downstream working section (m) 
    b  = 0.035;    % Radius of duct (m)
    a  = 0.0175;   % Radius of flame-holder (m)
    
    
    % Compute areas:
    A1 = pi*(b^2 - a^2);
    A2 = pi*b^2;
    
    % Spatial discretization:
    N = 50;
    [r,D] = fn_spatial_disc(a,b,N,order);
    
    % Forcing speakers?
    Geom.Speaker            =   false;
    Geom.Speaker_Downstream =   false;
    
    % Save to structure:
    Geom.Lu = Lu;
    Geom.Lb = Lb;
    Geom.a = a;
    Geom.b = b;
    
    Geom.A1 = A1;
    Geom.A2 = A2;
    
    Geom.D = D;
    Geom.r = r;
    Geom.N = N;
    
    
    Geom = orderfields(Geom);
end
