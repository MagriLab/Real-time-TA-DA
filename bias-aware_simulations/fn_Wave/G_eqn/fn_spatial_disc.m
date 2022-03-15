function [x,D] = fn_spatial_disc(a,b,N,order)
    % Spatial discretization:
    dx = (b-a)/(N-1);
    x = linspace(a,b,N)';
    % Diff matrix only for interior points:
    if order == 1
        D = zeros(N,N);
        for j = 1:N
            if j == 1
                D(j,j:j+1) = [1,-1]; % Forward difference
            else
                D(j,j-1:j) = [-1,1]; % Backward differentiation
            end
        end
        D = D/dx;
    elseif order == 2
        % Higher order:
        D = zeros(N,N);
        for j = 1:N
            if j == 1
                D(j,j:j+2) = [-3/2,2,-1/2]; % Forward difference
            elseif j == 2
                D(j,j-1:j+1) = [-1/2,0,1/2]; % central difference
            else
                D(j,j-2:j) = [1/2,-2,3/2]; % Backward differentiation
            end
        end
        D = D/dx;
    end
end