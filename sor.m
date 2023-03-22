function [phi] = sor(Nx, Ny, Ac, Ae, Aw, An, As, Bc, phi, n_sor)
for i_sor = 1:n_sor
    residual = 0;
    for i = 2:Ny + 1
        for j = 2:Nx + 1
            omega = (Ae(i, j)*phi(i, j + 1) + Aw(i, j)*phi(i, j - 1) + An(i, j)*phi(i + 1, j) + As(i, j)*phi(i - 1, j) + Bc(i,j))/Ac(i, j);
            residual = residual + abs(phi(i, j) - omega);
            phi(i, j) = omega;            
        end
    end
    residual = residual/(Nx*Ny);
    if residual < 1.0e-6
        break
    end
end
end

