nl = 5;
N = 20;
R = 1.5;
sz = 1;
I = 300;
mo = 4*pi*1e-7;
km = mo*I/(4*pi);
rw = 0.2;

% Angulo para cada punto de la espira
dtheta = 2*pi / N;
ang = 0:dtheta:(2*pi - dtheta);

% Inicializacion de posiciones e incrementos
s = 1;
for i = 1:nl
    Px(s:s+N-1) = R * cos(ang);
    Py(s:s+N-1) = R * sin(ang);
    Pz(s:s+N-1) = -nl/2*sz + (i-1)*sz;

    dx(s:s+N-1) = -Py(s:s+N-1) * dtheta;
    dy(s:s+N-1) =  Px(s:s+N-1) * dtheta;

    s = s + N;
end
dz = zeros(1, N*nl);

% Visualizacion de la espira
figure(1);
quiver3(Px, Py, Pz, dx, dy, dz, 0.5, '-r', 'LineWidth', 2);
view(-34, 33);
title('Corriente en espiras');
xlabel('x'); ylabel('y'); zlabel('z');
axis equal;

% PARTE 2: CAMPO MAGNETICO

% Definicion del espacio de calculo
ds = 0.1;
x = -5:ds:5; y = x; z = x;
Lx = length(x); Ly = length(y); Lz = length(z);

% Incializacion de los componentes del campo
dBx = zeros(Lx, Ly, Lz);
dBy = zeros(Lx, Ly, Lz);
dBz = zeros(Lx, Ly, Lz);

% Calculo del campo en cada punto del espacio 
for i= 1:Lx
    for j= 1:Ly
        for k= 1:Lz
            for l= 1:nl*N
                rx = x(i) - Px(l);
                ry = y(j) - Py(l);
                rz = z(k) - Pz(l);

                r = sqrt(rx^2 + ry^2 + rz^2 + rw^2); %Se suma rw^2 para evitar singularidades
                r3 = r^3;

                %Ley de Biot-Savart
                dBx (i,j,k) = dBx(i, j, k) + km * dy(l) * rz / r3;
                dBy (i,j,k) = dBy(i, j, k) + km * dx(l) * rz / r3;
                dBz (i,j,k) = dBz(i, j, k) + km * (dx(l) *ry - dy(l) * rx) / r3;
            end 
        end 
    end 
end

%Magnitud del campo 
Bmag= sqrt(dBx.^2 + dBy.^2 + dBz.^2);

%Corte en el plano XZ (y=0)
centery = round(Ly/2);
Bx_xz = squeeze(dBx(:, centery, :));
Bz_xz = squeeze(dBz(:, centery, :));
Bxz = squeeze(Bmag(:, centery, :));

%Visualización del campo 
figure(2)
hold on
pcolor (x,z, (Bxz').^ (1/3)); shading interp; colormap jet; colorbar 
h1 = streamslice(x, z, Bx_xz', Bz_xz', 3);
set(h1, "Color", [0.8, 1, 0.9]);
xlabel("x"); ylabel("z");
title ("Campo magnético generado por un solenoide");