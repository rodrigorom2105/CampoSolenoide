function [Bz, z] = campoB(ds, km, Px, Py, Pz, dx, dy, nl, N, rw, plot_option)
% Esta función calcula el campo magnético generado por un conjunto de espiras
% usando la Ley de Biot-Savart. Si plot_option es 1, grafica el campo en el plano XZ.
% Si es 0, devuelve el perfil del campo axial Bz a lo largo de z.

    % Rango de evaluación en el eje z con paso ds
    z = -5.2 : ds : 5.2;

    % Si se desea graficar, se usa una cuadrícula simétrica en x, y, z
    if plot_option
        x = z;
        y = z;
    else
        % Si no se grafica, se usa una pequeña región centrada en el eje z
        x = -0.1 : 0.01 : 0.1;
        y = -0.1 : 0.01 : 0.1;
    end

    % Tamaño del espacio tridimensional de cálculo
    Lx = length(x); Ly = length(y); Lz = length(z);

    % Inicialización de matrices para los componentes del campo magnético
    dBx = zeros(Lx, Ly, Lz, "single");
    dBy = zeros(Lx, Ly, Lz, "single");
    dBz = zeros(Lx, Ly, Lz, "single");

    % Cálculo del campo magnético en cada punto del espacio
    for i = 1:Lx
        for j = 1:Ly
            for k = 1:Lz
                for l = 1:nl*N
                    % Vector desde el elemento de corriente (Px,Py,Pz) al punto de evaluación (x,y,z)
                    rx = x(i) - Px(l);
                    ry = y(j) - Py(l);
                    rz = z(k) - Pz(l);

                    % Distancia desde el punto de corriente al punto de evaluación
                    % Se suma rw^2 para evitar divisiones por cero (singularidad)
                    r = sqrt(rx^2 + ry^2 + rz^2 + rw^2);
                    r3 = r^3;

                    % Aplicación de la Ley de Biot-Savart para cada componente del campo
                    dBx(i,j,k) = dBx(i,j,k) + km * dy(l) * rz / r3;
                    dBy(i,j,k) = dBy(i,j,k) + km * dx(l) * rz / r3;
                    dBz(i,j,k) = dBz(i,j,k) + km * (dx(l) * ry - dy(l) * rx) / r3;
                end
            end
        end
    end

    % Si se desea graficar el campo en 2D
    if plot_option
        % Magnitud total del campo magnético
        Bmag = sqrt(dBx.^2 + dBy.^2 + dBz.^2);

        % Se toma un corte en el plano Y=0 (centro del espacio en Y)
        centery = round(Ly / 2);
        Bx_xz = squeeze(dBx(:, centery, :));   % componente x del campo en plano XZ
        Bz_xz = squeeze(dBz(:, centery, :));   % componente z del campo en plano XZ
        Bxz = squeeze(Bmag(:, centery, :));    % magnitud del campo en plano XZ

        % Visualización del campo magnético
        figure
        hold on
        % Mapa de colores del campo magnético, se le aplica raíz cúbica para mejor contraste visual
        pcolor(x, z, (Bxz').^(1/3));
        shading interp;
        colormap jet;
        colorbar;

        % Añade líneas de flujo del campo (streamlines)
        h1 = streamslice(x, z, Bx_xz', Bz_xz', 1);
        set(h1, "Color", [0.8, 1, 0.9]);  % Color verde claro

        % Etiquetas y título del gráfico
        xlabel("x"); ylabel("z");
        title("Campo magnético generado por un solenoide");
    else
        % Si no se grafica, se calcula el perfil axial del campo Bz a lo largo de z
        idx_x = ceil(Lx / 2);
        idx_y = ceil(Ly / 2);

        % Se toma el valor de Bz en el centro (x=0, y=0) para todas las z
        Bz = squeeze(dBz(idx_x, idx_y, :));

        % Se calcula el gradiente del campo Bz respecto a z
        dBz_dz_profile = diff(Bz) ./ diff(z);
        z_mid = z(1:end-1) + diff(z)/2;  % puntos medios de z para el gradiente

        % Gráfica del gradiente de Bz
        figure
        grid on;
        plot(z_mid, dBz_dz_profile, "r-", "LineWidth", 2);
        xlabel('z'); ylabel('Valor');
        title('Gradiente del campo');
    end

end
