function postproc_centerline(Nx, Ny, y, u)
figure;
plot(u(2:Ny + 1, round(Nx/2)), y);
xlabel('u');
ylabel('x');
title('Centerline u-velocity');
end

