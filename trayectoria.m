function trayectoria(Bz, z, mag, m, zo, dt, vz, gamma)
% Esta función simula y grafica la trayectoria de un imán que cae bajo la gravedad
% y la fuerza magnética inducida por un campo Bz, comparándola con la caída libre.
%
% Parámetros:
% Bz     -> Perfil del campo magnético (vector)
% z      -> Coordenadas z correspondientes a Bz
% mag    -> Momento magnético del imán
% m      -> Masa del imán
% zo     -> Posición inicial en z
% dt     -> Intervalo de tiempo
% vz     -> Velocidad inicial
% gamma  -> Coeficiente de fricción viscosa

    % Fuerza peso (constante): w = -mg
    w = -m * 9.81;

    % Inicialización de vectores de posición y velocidad
    zm(1) = zo;          % Posición del imán con campo magnético
    zmfree(1) = zo;      % Posición del imán sin campo (caída libre)
    vz(1) = 0.7;         % Velocidad inicial con campo
    vzfree(1) = 0;       % Velocidad inicial sin campo
    tt(1) = 0;           % Tiempo inicial
    cc = 1;              % Contador de pasos

    % Bucle de simulación: continúa hasta que el imán cae por debajo de -3 metros
    while zm(cc) > -3
        % Cálculo numérico del gradiente de Bz usando diferencias centrales
        delta = 0.005;
        Bz_forward  = interp1(z, Bz, zm(cc) + delta, "linear", "extrap");
        Bz_backward = interp1(z, Bz, zm(cc) - delta, "linear", "extrap");
        dBz_dz = (Bz_forward - Bz_backward) / (2 * delta);

        % Cálculo de la fuerza magnética: Fm = -μ * dBz/dz
        Fm(cc) = -mag * dBz_dz;

        % Fuerza de fricción viscosa: proporcional a la velocidad
        Ff = -gamma * vz(cc);

        % Fuerza total: magnética + peso + fricción
        F(cc) = Fm(cc) + w + Ff;

        % Aceleración instantánea: F = m*a => a = F/m
        a = F(cc) / m;

        % Movimiento libre (sin campo)
        zmfree(cc+1) = zmfree(cc) + vzfree(cc)*dt + 0.5*(-9.81)*dt^2;
        vzfree(cc+1) = vzfree(cc) - 9.81 * dt;

        % Movimiento con campo magnético
        zm(cc+1) = zm(cc) + vz(cc)*dt + 0.5*a*dt^2;
        vz(cc+1) = vz(cc) + a*dt;

        % Actualizar tiempo
        tt(cc+1) = tt(cc) + dt;

        % Avanzar contador
        cc = cc + 1;

        % Condición de paro: si la velocidad es muy pequeña (reposo virtual)
        if abs(vz(cc)) < 1e-3
            break;
        end
    end

    % Visualización: comparación entre caída libre y caída con campo magnético
    figure;
    plot(tt, zm, '-r', 'LineWidth', 2); hold on;               % Con campo
    plot(tt, zmfree, '--b', 'LineWidth', 2);                   % Sin campo
    legend('Fall over a current ring', ...
           'Free fall (no Magnetic force)', ...
           'Location', 'best');
    xlabel('time (s)');
    ylabel('z position (m)');
    title('Position vs time of a Magnetic dipole falling through a current ring');
    grid on;
    ylim([-6 6]);  % Ajuste del eje y para mantener escala fija
end
