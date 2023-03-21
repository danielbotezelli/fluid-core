function res = residual(Nx, Ny, p_prime, Ac_u, dy)
res = 0;
for i = 2:Ny + 1
    for j = 2:Nx + 1
        omega = abs(0.5*(p_prime(i, j - 1) - p_prime(i, j + 1))/Ac_u(i, j)*dy);
        res = res + omega;          
    end
end
res = res/(Nx*Ny);
end

