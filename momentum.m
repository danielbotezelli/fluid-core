function [Ac_u, Ac_v, Ae, Aw, An, As, Bc_u, Bc_v] = momentum(Nx, Ny, dx, dy, Ac_u, Ac_v, Ae, Aw, An, As, Bc_u, Bc_v, ...
                                                     u, v, uf, vf, grad_p_i, grad_p_j, rho, mu, lambda, U, channelflow_model, liddriven_model)
% Fully developed gradient
if channelflow_model
    grad_ub_i = 0;
    grad_ub_j = 0;
    uf(:, Nx + 1) = u(:, Nx + 1) + 0.5*(grad_ub_i*dx + grad_ub_j*dy);
    u(:, Nx + 2) = u(:, Nx + 1);
    uf(:, Nx) = u(:, Nx + 1);
    
    v(:, Nx + 2) = v(:, Nx + 1);
    vf(1:Ny + 1, Nx + 2) = v(1:Ny + 1, Nx + 1);
end
for i = 2:Ny + 1
    for j = 2:Nx + 1
        %% Matrix
        % Mass flow
        me = rho*uf(i, j)*dy;
        mw = rho*uf(i, j - 1)*dy;
        mn = rho*vf(i, j)*dx;
        ms = rho*vf(i - 1, j)*dx;

        % Linkers
        if j < Nx + 1,  Ae(i, j) = max(0.0, -me) + mu*dy/dx;
        else,           Ae(i, j) = 0;
        end

        if j > 2,       Aw(i, j) = max(0.0, mw) + mu*dy/dx;
        else,           Aw(i, j) = 0;
        end
        
        if i < Ny + 1,  An(i, j) = max(0.0, -mn) + mu*dx/dy;
        else,           An(i, j) = 0;
        end

        if i > 2,       As(i, j) = max(0.0, ms) + mu*dx/dy;
        else,           As(i, j) = 0;
        end

        Ac_u(i, j) = Ae(i, j) + Aw(i, j) + An(i, j) + As(i, j) + (me - mw) + (mn - ms);
        Ac_v(i, j) = Ae(i, j) + Aw(i, j) + An(i, j) + As(i, j) + (me - mw) + (mn - ms);

        % Face gradient
        grad_uf_i_e = (u(i, j + 1) - u(i, j))/dx;
        grad_uf_i_w = -(u(i, j - 1) - u(i, j))/dx;
        grad_uf_i_n = 0;
        grad_uf_i_s = 0;
        
        grad_uf_j_e = 0;
        grad_uf_j_w = 0;
        grad_uf_j_n = (u(i + 1, j) - u(i, j))/dy;
        grad_uf_j_s = -(u(i - 1, j) - u(i, j))/dy;

        grad_vf_i_e = (v(i, j + 1) - v(i, j))/dx;
        grad_vf_i_w = -(v(i, j - 1) - v(i, j))/dx;
        grad_vf_i_n = 0;
        grad_vf_i_s = 0;
        
        grad_vf_j_e = 0;
        grad_vf_j_w = 0;
        grad_vf_j_n = (v(i + 1, j) - v(i, j))/dy;
        grad_vf_j_s = -(v(i - 1, j) - v(i, j))/dy;

        if j == Nx + 1, grad_uf_i_e = 2*grad_uf_i_e; end
        if j == 2,      grad_uf_i_w = 2*grad_uf_i_w; end
        if i == Ny + 1, grad_uf_j_n = 2*grad_uf_j_n; end
        if i == 2,      grad_uf_j_s = 2*grad_uf_j_s; end

        % High resolution
        ue_HR = 0.5*(grad_uf_i_e*dx + grad_uf_j_e*dy);
        uw_HR = -0.5*(grad_uf_i_w*dx + grad_uf_j_w*dy);
        un_HR = 0.5*(grad_uf_i_n*dx + grad_uf_j_n*dy);
        us_HR = -0.5*(grad_uf_i_s*dx + grad_uf_j_s*dy);

        ve_HR = 0.5*(grad_vf_i_e*dx + grad_vf_j_e*dy);
        vw_HR = -0.5*(grad_vf_i_w*dx + grad_vf_j_w*dy);
        vn_HR = 0.5*(grad_vf_i_n*dx + grad_vf_j_n*dy);
        vs_HR = -0.5*(grad_vf_i_s*dx + grad_vf_j_s*dy);
        
        % Source
        Bc_u(i, j) = -(me*ue_HR - mw*uw_HR + mn*un_HR - ms*us_HR) ...
                   + mu*(grad_uf_i_e*dy - grad_uf_i_w*dy + grad_uf_j_n*dx - grad_uf_j_s*dx) ...
                   - grad_p_i(i, j)*dx*dy;

        Bc_v(i, j) = -(me*ve_HR - mw*vw_HR + mn*vn_HR - ms*vs_HR) ...
                   + mu*(grad_vf_i_e*dy - grad_vf_i_w*dy + grad_vf_j_n*dx - grad_vf_j_s*dx) ...
                   - grad_p_j(i, j)*dx*dy;
        
        %% Boundary conditions
        % East
        if j == Nx + 1
            m_dot = 0;
            if channelflow_model
                m_dot = rho*uf(i, j)*dy;
            end
            f_tau = 0;
            if liddriven_model  
                f_tau = 2*mu*dx/dy;
            end
            Ac_u(i, j) = Ac_u(i, j) + m_dot + 2*mu*dy/dx;
            Ac_v(i, j) = Ac_v(i, j) + m_dot + f_tau + 2*mu*dy/dx;
            Bc_u(i, j) = Bc_u(i, j) + m_dot*uf(i, j) + 2*mu*dy/dx*uf(i, j);
            Bc_v(i, j) = Bc_v(i, j) + m_dot*vf(i, j) + 2*mu*dy/dx*vf(i, j);
        end

        % West
        if j == 2
            af = 0;
            if channelflow_model
                af = max(0.0, mw) + 2*mu*dx/dy;
            end
            f_tau = 0;
            if liddriven_model  
                f_tau = 2*mu*dx/dy;
                af = max(0.0, mw) + 2*mu*dx/dy;
            end
            Ac_u(i, j) = Ac_u(i, j) + af;
            Ac_v(i, j) = Ac_v(i, j) + f_tau + af;
            Bc_u(i, j) = Bc_u(i, j) + af*U;
            Bc_v(i, j) = Bc_v(i, j);
        end

        % North
        if i == Ny + 1
            f_tau = 0;
            if or(channelflow_model, liddriven_model)
                f_tau = 2*mu*dy/dx;
            end
            momentum = 0;
            if liddriven_model
                momentum = f_tau*U + 2*mu*dy/dx;
            end
            Ac_u(i, j) = Ac_u(i, j) + f_tau + 2*mu*dy/dx;
            Ac_v(i, j) = Ac_v(i, j) + 2*mu*dy/dx;
            Bc_u(i, j) = Bc_u(i, j) + momentum;
            Bc_v(i, j) = Bc_v(i, j);     
        end

        % South
        if i == 2
            f_tau = 0;
            if or(channelflow_model, liddriven_model)  
                f_tau = 2*mu*dy/dx;
            end
            Ac_u(i, j) = Ac_u(i, j) + f_tau + 2*mu*dy/dx;
            Ac_v(i, j) = Ac_v(i, j) + 2*mu*dy/dx;
            Bc_u(i, j) = Bc_u(i, j);
            Bc_v(i, j) = Bc_v(i, j);
        end
        
        %% Relaxation
        Ac_u(i, j) = Ac_u(i, j)/lambda;
        Ac_v(i, j) = Ac_v(i, j)/lambda;
        Bc_u(i, j) = Bc_u(i, j) + (1 - lambda)*Ac_u(i, j)*u(i, j);
        Bc_v(i, j) = Bc_v(i, j) + (1 - lambda)*Ac_v(i, j)*v(i, j);
    end
end
end

