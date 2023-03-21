function [grad_phi_i, grad_phi_j] = grad(Nx, Ny, phi, dx, dy, grad_phi_i, grad_phi_j, weight)
for i = 2:Ny + 1
    for j = 2:Nx + 1
        wx = 1/((dx^2)^(weight/2));
        wy = 1/((dy^2)^(weight/2));
%         w = 1;
        a0 = 2*wx*dx^2;
        a1 = 0;
        a2 = 0;
        a3 = 2*wy*dy^2;
        b0 = wx*dx*(phi(i, j + 1) - phi(i, j - 1));
        b1 = wy*dy*(phi(i + 1, j) - phi(i - 1, j));

        % Determinant
        det = a0*a3 - a1*a2;
        det_x = b0*a3 - a1*b1;
        det_y = a0*b1 - b0*a2;

        % Gradient
        if (det ~= 0)
            grad_phi_i(i, j) = det_x/det;
            grad_phi_j(i, j) = det_y/det;
        else
            grad_phi_i(i, j) = 0;
            grad_phi_j(i, j) = 0;
        end
    end
end
end

