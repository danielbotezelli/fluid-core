function [Ac_p, Ae, Aw, An, As, Bc, p_prime] = continuity(Nx, Ny, dx, dy, Ac_u, Ac_v, Ac_p, Ae, Aw, An, As, Bc, p_prime, uf, vf, rho, liddriven_model)
for i = 2:Ny + 1
    for j = 2:Nx + 1
        Df_u_e = 0.5*(dx*dy/Ac_u(i, j) + dx*dy/Ac_u(i, j + 1));
        Df_u_w = 0.5*(dx*dy/Ac_u(i, j) + dx*dy/Ac_u(i, j - 1));

        Df_v_n = 0.5*(dx*dy/Ac_v(i, j) + dx*dy/Ac_v(i + 1, j));
        Df_v_s = 0.5*(dx*dy/Ac_v(i, j) + dx*dy/Ac_v(i - 1, j));

        Df_e = (Df_u_e*dy)^2/(dx*Df_u_e*dy);
        Df_w = (Df_u_w*dy)^2/(dx*Df_u_w*dy);
        Df_n = (Df_v_n*dx)^2/(dy*Df_v_n*dx);
        Df_s = (Df_v_s*dx)^2/(dy*Df_v_s*dx);

        if j < Nx + 1,  Ae(i, j) = rho*Df_e;
        else,           Ae(i, j) = 0; 
        end

        if j > 2,       Aw(i, j) = rho*Df_w;
        else,           Aw(i, j) = 0; 
        end

        if i < Ny + 1,  An(i, j) = rho*Df_n;
        else,           An(i, j) = 0; 
        end

        if i > 2,       As(i, j) = rho*Df_s;
        else,           As(i, j) = 0; 
        end

        Ac_p(i, j) = Ae(i, j) + Aw(i, j) + An(i, j) + As(i, j);

        % Source
        Bc(i, j) = - rho*(uf(i, j) - uf(i, j - 1))*dy - rho*(vf(i, j) - vf(i - 1, j))*dx;

        % Redefine
        p_prime(i, j) = 0;
        
        % Anchor
        if liddriven_model
            Ae(2, 2) = 0;
            Aw(2, 2) = 0;
            An(2, 2) = 0;
            As(2, 2) = 0;
            Ac_p(2, 2) = 1;
            Bc(2, 2) = 0;
        end
    end
end
end

