function [ D, g ] = Cheb(N)
%   Function that computes the Chebyshev collocation derivative matrix (D) 
%   and the Chevyshev grid of (N + 1) points.  

    g   =   cos(pi * (0:N)/N)';
    c   =   [2; ones(N-1,1); 2].*(-1).^(0:N)';
    X   =   repmat(g,1,N+1);
    dX  =   X - X';
    
    D   =   (c * (1./ c)')./ (dX + (eye(N+1)));
    D   =   D - diag(sum(D'));  
end

