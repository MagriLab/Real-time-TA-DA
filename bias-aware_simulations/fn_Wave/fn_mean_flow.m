function [Mean] = fn_mean_flow(u1,T1,p2,Qbar,Geom, R_in, R_out)
    % Function that compute the Mean flow for the entire geometry. 
    
    %  Ideal gas properties for air:
    gamma = 1.4;
    R_gas = 287.1;
    Cp = 1005;
    
    % Geometric definitions:
    A1 = Geom.A1;
    A2 = Geom.A2;
    a  = Geom.a;
    b  = Geom.b;
    r  = Geom.r;
   
    
    
%     if u1 > 1
        % Solve for the mean flow:
        x_init = [p2,1,T1+Qbar/(A1*u1*Cp),1,u1*A1/A2];
        [xs,~,EF] = fsolve(@(x) Flux_Balance(x,A2,A1,Cp,R_gas,u1,T1,p2,Qbar),...
            x_init,optimoptions('fsolve','Display','none'));
        if EF > 0
           p1 = xs(1);  rho1 = xs(2);  T2 = xs(3); rho2 = xs(4);  u2 = xs(5); 
        else
            error('base flow error')
        end
        % Compute properties before the jump
        c1 = sqrt(gamma*R_gas*T1);
        M1 = u1/c1;
        % Compute properties after the jump
        c2 = sqrt(gamma*R_gas*T2);
        M2 = u2/c2;
        % Time delays
        Lu = Geom.Lu;
        Lb = Geom.Lb;
        Tu = 2*Lu/(c1 * (1 - M1^2));
        Td = 2*Lb/(c2 * (1 - M2^2));
        Tau  = 0.4*Lb/u1;                       % Time delay of the flame
%     else
%         % Zero mach number limit
%         p1      =   p2;
%         T2      =   T1;
%         u2      =   u1 * A1/A2;
%         rho1    =   p1 / (R_gas * T1);
%         rho2    =   p2 / (R_gas * T2);
%         % Compute properties before the jump
%         c1 = sqrt(gamma*R_gas*T1);
%         M1 = u1/c1;
%         % Compute properties after the jump
%         c2 = sqrt(gamma*R_gas*T2);
%         M2 = u2/c2;
%         % Time delays
%         Tu  =   0.0069;
%         Td  =   0.0034;
%         Tau =   0.01;
%     end
    
    if contains(Geom.BC1,'closed')
        R_in	=   (1-M1)/(1+M1);
    end
    
    % Flame parameters:
    Su   = 0.09*u1;                         % Flame speed
    Abar = pi*(b^2 - a^2) * u1/Su;          % Mean flame area
    xibar = (r-a)*(u1^2 - Su^2)^(1/2)/Su;   % Mean flame front
    
    
    % X, Y Matrices:
%     X(1,1) = (M2*c2*(M1 - 1))/c1 - M1*(M1 - 2) - A2/A1; 
%     X(1,2) = (A2*(M2 + 1))/A1; 
%     X(2,1) = (M2^2*c2^2*(M1 - 1))/(2*c1^2) - (M1*(M1 - 1)*(M1 - 2))/2 - (M1 - 1)/(gamma - 1); 
%     X(2,2) = (A2*c2*(M2 + 1)*(M2*gamma - M2 + 1))/(A1*c1*(gamma - 1)); 
%     Y(1,1) = (M2*c2*(M1 - 1))/c1 - ((M1 - 1)*(A2 + A1*M1^2 + 2*A1*M1))/(A1*(M1 + 1)); 
%     Y(1,2) = (A2*R_out*(M2 - 1))/A1; 
%     Y(2,1) = (M2^2*c2^2*(M1 - 1))/(2*c1^2) - (M1*(M1 - 1)*(M1 + 2))/2 - (M1 - 1)/(gamma - 1); 
%     Y(2,2) = -(A2*R_out*c2*(M2 - 1)*(M2 - M2*gamma + 1))/(A1*c1*(gamma - 1));  
    
    X(1,1) = (M2*c2*(M1 - 1))/c1 - M1*(M1 - 2) - A2/A1; 
    X(1,2) = (A2*(M2 + 1))/A1; 
    X(2,1) = (M2^2*c2^2*(M1 - 1))/(2*c1^2) - (M1*(M1 - 1)*(M1 - 2))/2 - (M1 - 1)/(gamma - 1); 
    X(2,2) = (A2*c2*(M2 + 1)*(M2*gamma - M2 + 1))/(A1*c1*(gamma - 1)); 
    Y(1,1) = (A2*R_in)/A1 + M1*R_in*(M1 + 2) - (M2*R_in*c2*(M1 + 1))/c1; 
    Y(1,2) = (A2*R_out*(M2 - 1))/A1; 
    Y(2,1) = (R_in*(M1 + 1))/(gamma - 1) + (R_in*(M1 + 1)*(2*M1*c1^2 + M1^2*c1^2 - M2^2*c2^2))/(2*c1^2); 
    Y(2,2) = -(A2*R_out*c2*(M2 - 1)*(M2 - M2*gamma + 1))/(A1*c1*(gamma - 1)); 
    
    % Speaker Forcing:
    F(1,1) = (c1 + M1*c1 - M2*c2)/A1 - (c1*(A1 - A2))/(A1^2*(M1 + 1));
    F(2,1) = c1/(A1*(gamma - 1)) - (M2^2*c2^2)/(2*A1*c1) + (M1*c1*(M1 + 2))/(2*A1);
    % Downstream speaker forcing:
    D(1,1) = (A2*(M2 - 1))/A1;
    D(2,1) = -(A2*c2*(M2 - 1)*(M2 - M2*gamma + 1))/(A1*c1*(gamma - 1));
    
    % Save to structure:
    Mean.M1 = M1; 
    Mean.T1 = T1; 
    Mean.c1 = c1;
    Mean.u1 = u1;
    Mean.p1 = p1;
    Mean.rho1 = rho1;
    
    Mean.M2 = M2; 
    Mean.T2 = T2; 
    Mean.c2 = c2;
    Mean.u2 = u2;
    Mean.p2 = p2;
    Mean.rho2 = rho2;
    
    Mean.Qbar = Qbar;
    
    Mean.gamma = gamma;
    Mean.Cp = Cp;
    Mean.R_gas = R_gas;
    
    Mean.Tau = Tau;
    Mean.Su = Su;   
    Mean.Abar = Abar ;

    Mean.xibar = Unstable_Flame();
%     Mean.xibar = xibar;
    
    Mean.Tu = Tu;
    Mean.Td = Td;
    Mean.Su = Su;
    
    Mean.X = X;
    Mean.Y = Y;
    Mean.F = F;
    Mean.D = D;
    Mean.IX = inv(X);
    Mean.IXY = X\Y;
    
    Mean.R_in   =   R_in;
    Mean.R_out	=   R_out;
    
    Mean = orderfields(Mean);
end

%% Flux Equations
function [Eq] = Flux_Balance(x,A2,A1,Cp,R_gas,u1,T1,p2,Q)
    % Function to balance the fluxes and compute the mean flow on the next
    % side of the jump condition assuming:
    %          m2 - m1    = 0,
    %         fx2 - fx1   = (A2 - A1)*p1,
    %       m2*H2 - m1*H1 = Q,
    
    % Set the unknowns:
    p1 = x(1);  rho1 = x(2); T2 = x(3); rho2 = x(4); u2 = x(5);
    % Compute properties before the jump:
    m1  = A1*rho1*u1;
    fx1 = A1*p1 + m1*u1;
    H1  = Cp*T1 + 1/2*u1^2;
    % Compute properties of after the jump:
    m2  = A2*rho2*u2;
    fx2 = A2*p2 + m2*u2;
    H2  = Cp*T2 + 1/2*u2^2;
    % Set the equations:
    % Ideal Gas:
    Eq(1) = 1 - (T1*rho1*R_gas)/p1;
    Eq(2) = 1 - (T2*rho2*R_gas)/p2;
    % Mass conservation:
    Eq(3) = 1 - m2/m1;
    % Momentum conservation:
    Eq(4) = 1 - fx2/(fx1 + (A2 - A1)*p1);
    % Energy conservation:
    Eq(5) = 1 - m2*H2/(m1*H1 + Q);
end

function [Xi_unst] = Unstable_Flame()
    Xi_unst = [0.00540448322988493
        0.0187036043186722
        0.0319627241581763
        0.0452562786905994
        0.0587329486729451
         0.072585690754649
        0.0867599638162835
         0.100971576709858
          0.11357010731261
         0.124179874053602
         0.129545706496416
         0.132953721400073
         0.134856306427151
         0.137165347927167
         0.139366697398413
         0.141020162024224
         0.141910076520878
         0.142662446269569
         0.142734499037492
         0.142810998064086
         0.142876837519442
          0.14294663411858
          0.14302008645732
         0.143097047414309
         0.143177593696378
         0.143261748589889
         0.143349533811358
          0.14344097470085
         0.143536097891067
         0.143634931294008
         0.143737504213545
         0.143843847410752
         0.143953993147733
          0.14406797524119
         0.144185829120302
         0.144307591888205
         0.144433302387561
         0.144563001270522
         0.144696731073433
         0.144834536296601
         0.144976463489527
         0.145122561341999
         0.145272880781523
         0.145427475077584
         0.145586399953283
         0.145749713704972
         0.145917477330535
         0.146089754667063
         0.146266612538734
         0.146448120915783];
end