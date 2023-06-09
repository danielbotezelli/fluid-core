function [u, v, uf, vf, p] = pressure_correction(Nx, Ny, u, v, uf, vf, p, grad_p_i, grad_p_j, p_prime, grad_p_prime_i, grad_p_prime_j,...
                                                 Ac_u, Ac_v, lambda_p, dx, dy, channelflow_model, rho, U, Ly)
for i = 2:Ny + 1
    for j = 2:Nx + 1
        % Pressure
        p(i, j) = p(i, j) + lambda_p*p_prime(i, j);

        % Velocity
        u(i, j) = u(i, j) - dx*dy/Ac_u(i, j)*grad_p_prime_i(i, j);
        v(i, j) = v(i, j) - dx*dy/Ac_v(i, j)*grad_p_prime_j(i, j);
    end
end

% Face u-velocity
for i = 2:Ny + 1
    for j = 2:Nx
        Df_u_e = 0.5*(dx*dy/Ac_u(i, j) + dx*dy/Ac_u(i, j + 1));
        grad_p_prime_i_e = (p_prime(i, j + 1) - p_prime(i, j))/dx;

        uf(i, j) = uf(i, j) - Df_u_e*grad_p_prime_i_e;
    end
end

% Face v-velocity
for i = 2:Ny
    for j = 2:Nx + 1
        Df_v_n = 0.5*(dx*dy/Ac_v(i, j) + dx*dy/Ac_v(i + 1, j));
        grad_p_prime_j_n = (p_prime(i + 1, j) - p_prime(i, j))/dy;

        vf(i, j) = vf(i, j) - Df_v_n*grad_p_prime_j_n;
    end
end

% Boundary pressure
p(2:Ny + 1, Nx + 2) = p(2:Ny + 1, Nx + 1) + 0.5*grad_p_i(2:Ny + 1, Nx + 1)*dx;  % East
p(2:Ny + 1, 1) = p(2:Ny + 1, 2) - 0.5*grad_p_i(2:Ny + 1, 2)*dx;                 % West 
p(Ny + 2, 2:Nx + 1) = p(Ny + 1, 2:Nx + 1) + 0.5*grad_p_j(Ny + 1, 2:Nx + 1)*dy;  % North
p(1, 2:Nx + 1) = p(2, 2:Nx + 1) - 0.5*grad_p_j(2, 2:Nx + 1)*dy;                 % South

% Outlet
if channelflow_model
    mass_in = rho*U*Ly;
    mass_out = 0;
    for i = 1:Ny + 2
        mass_out = mass_out + rho*u(i, Nx + 2)*dy;
    end
    if mass_out ~= 0
        u(:, Nx + 2) = u(:, Nx + 2)*(mass_in/mass_out);
        v(:, Nx + 2) = v(:, Nx + 2)*(mass_in/mass_out);
        
        uf(:, Nx + 1) = u(:, Nx + 2);
        vf(1:Ny + 1, Nx + 2) = v(1:Ny + 1, Nx + 2);
    end
end
end

