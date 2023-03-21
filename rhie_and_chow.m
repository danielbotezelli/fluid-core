function [u, v, uf, vf] = rhie_and_chow(Nx, Ny, u, v, uf, vf, p, grad_p_i, grad_p_j, ...
                                        Ac_u, Ac_v, dx, dy, channelflow_model)
for i = 2:Ny + 1
    for j = 2:Nx
        grad_pf_bar_e = 0.5*(grad_p_i(i, j) + grad_p_i(i, j + 1));
        grad_pf_e = (p(i, j + 1) - p(i, j))/dx - grad_pf_bar_e;
        Df_e = 0.5*(dx*dy/Ac_u(i, j) + dx*dy/Ac_u(i, j + 1));

        uf(i, j) = 0.5*(u(i, j) + u(i, j + 1)) - Df_e*grad_pf_e;
    end
end

for i = 2:Ny
    for j = 2:Nx + 1
        grad_pf_bar_n = 0.5*(grad_p_j(i, j) + grad_p_j(i + 1, j));
        grad_pf_n = (p(i + 1, j) - p(i, j))/dy - grad_pf_bar_n;
        Df_n = 0.5*(dx*dy/Ac_v(i, j) + dx*dy/Ac_v(i + 1, j));

        vf(i, j) = 0.5*(v(i, j) + v(i + 1, j)) - Df_n*grad_pf_n;
    end
end

% Outlet
if channelflow_model  
    u(:, Nx + 2) = u(:, Nx + 1);
    uf(:, Nx + 1) = u(:, Nx + 1);
    uf(:, Nx) = u(:, Nx + 1);

    v(:, Nx + 2) = v(:, Nx + 1);
    vf(1:Ny + 1, Nx + 2) = v(1:Ny + 1, Nx + 1);
end
end
