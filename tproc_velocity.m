function tproc_velocity(u, v, X, Y, x, y, Nx, Ny, Lx, channelflow_model)
mycolormap = customcolormap(linspace(0,1,11), {'#68011d','#b5172f','#d75f4e','#f7a580','#fedbc9','#f5f9f3','#d5e2f0','#93c5dc','#4295c1','#2265ad','#062e61'});
contourf(X, Y, sqrt(u(2:Ny+1, 2:Nx+1).^2 + v(2:Ny+1, 2:Nx+1).^2), 10, 'LineStyle', 'none');
colormap(mycolormap);
colorbar;
hold on;
sc = round(Ny/15);
quiver(X(1:sc:end, 1:sc:end), Y(1:sc:end, 1:sc:end), u(2:sc:Ny + 1, 2:sc:Nx + 1), v(2:sc:Ny + 1, 2:sc:Nx + 1), 'Color', 'w', 'LineWidth', 0.7, 'AutoScale', 'on', 'AutoScaleFactor', 0.5);
if channelflow_model
    for i = 1:5
        hold on;
        plot(u(2:Ny + 1, round(i*Nx/5)) + (i*Lx - Lx)/5, y, 'Color', 'w', 'LineWidth', 3.0);
    end
end
xlabel('x');
ylabel('y');
title('Quiver');
daspect([1 1 1]);
axis equal;
box on;
set(gca, 'FontSize', 14);
set(gca, 'LineWidth', 2);
set(gca, 'Layer', 'top');
set(gca, 'Color', 'w');
set(gcf, 'Color', 'w');
set(gcf, 'WindowState', 'maximized');
xlim([x(1) x(end)]);
ylim([y(1) y(end)]);
drawnow;
end

