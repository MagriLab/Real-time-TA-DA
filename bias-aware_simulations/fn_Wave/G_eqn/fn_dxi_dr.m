function [dxi_dr] = fn_dxi_dr(xi,Ut,Mean,Geom)
    % Compute dxi/dr and applies BCs. 
    D = Geom.D;
    Su = Mean.Su;
    % Set the boundary condition:
    if Ut >= Su && xi(1) > -1e-10
        Gr = sqrt(Ut^2/Su^2 - 1);
    elseif Ut < Su || xi(1) < 0
        Gr = 0;
    else
         error('error')
    end
    dxi_dr = D*xi;
    dxi_dr(1) = Gr;
    
    if any(dxi_dr > 1e10)
        a=0
    end
end