function [grad_phi_i, grad_phi_j] = gradgauss_prime(Nx, Ny, phi, dx, dy, grad_phi_i, grad_phi_j)
for i = 2:Ny + 1
    for j = 2:Nx + 1
        phif_e = 0.5*(phi(i, j) + phi(i, j + 1));
        phif_w = 0.5*(phi(i, j) + phi(i, j - 1));
        phif_n = 0.5*(phi(i, j) + phi(i + 1, j));
        phif_s = 0.5*(phi(i, j) + phi(i - 1, j));
        
        if j == Nx + 1, phif_e = phi(i, j); end
        if j == 2,      phif_w = phi(i, j); end
        if i == Ny + 1, phif_n = phi(i, j); end
        if i == 2,      phif_s = phi(i, j); end

        grad_phi_i(i, j) = (phif_e*dy - phif_w*dy)/(dx*dy);
        grad_phi_j(i, j) = (phif_n*dx - phif_s*dx)/(dx*dy);       
    end
end
end