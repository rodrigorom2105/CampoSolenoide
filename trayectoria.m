function zm = trayectoria(Bz, z, mag, m, zo, dt, vz, gamma)
% Simula la caída de un imán a través de un campo magnético Bz usando el método de Runge-Kutta de orden 4.
% También compara con la caída libre (sin campo magnético).
%
% Parámetros:
% Bz     -> Perfil del campo magnético en z
% z      -> Vector de posiciones z correspondientes a Bz
% mag    -> Momento magnético del imán
% m      -> Masa del imán
% zo     -> Posición inicial
% dt     -> Paso de tiempo
% vz     -> Velocidad inicial con campo (aunque será sobrescrita)
% gamma  -> Coeficiente de fricción viscosa

    w = -m * 9.81;  % Fuerza peso (constante)

    % Inicialización de vectores para trayectoria
    zm(1) = zo;         % Posición con campo magnético
    zmfree(1) = zo;     % Posición en caída libre
    vz(1) = 0.7;        % Velocidad inicial (con campo)
    vzfree(1) = 0;      % Velocidad inicial en caída libre
    tt(1) = 0;          % Tiempo
    cc = 1;             % Contador de iteraciones

    % Simulación paso a paso
    while zm(cc) > -3
        % Movimiento libre (sin campo)
        zmfree(cc+1) = zmfree(cc) + vzfree(cc)*dt + 0.5*(-9.81)*dt^2;
        vzfree(cc+1) = vzfree(cc) - 9.81 * dt;

        % Movimiento con frenado y campo magnético usando Runge-Kutta 4
        z_axis = z;  % Copia del eje z para pasar a la función de aceleración
        [zm(cc+1), vz(cc+1)] = rk4_step(zm(cc), vz(cc), dt, @(z, v) a_total(z, v, Bz, z_axis, mag, gamma, m));

        % Actualización del tiempo
        tt(cc+1) = tt(cc) + dt;

        % Avanza el contador
        cc = cc + 1;

        % Condición de parada: velocidad muy baja (reposo)
        if abs(vz(cc)) < 1e-3
            break;  
        end
    end

    % Gráfica comparativa de trayectorias
    figure;
    plot(tt, zm, '-r', 'LineWidth', 2); hold on;
    plot(tt, zmfree, '--b', 'LineWidth', 2);
    legend('Fall over a current ring', 'Free fall (no Magnetic force)', 'Location', 'best');
    xlabel('time (s)');
    ylabel('z position (m)');
    title('Position vs time of a Magnetic dipole falling through a current ring');
    grid on;
    ylim([-6 6]);  % Limita el eje y para mantener escala fija

    % FUNCIÓN ANIDADA: Calcula la aceleración total del imán en función de su posición y velocidad
    function a = a_total(z, v, Bz, z_axis, mag, gamma, m)
        delta = 0.005;

        % Interpolación de Bz en puntos adyacentes para calcular derivada
        Bz_forward  = interp1(z_axis, Bz, z + delta, "linear", "extrap");
        Bz_backward = interp1(z_axis, Bz, z - delta, "linear", "extrap");
        dBz_dz = (Bz_forward - Bz_backward) / (2 * delta);  % Derivada central

        % Fuerza magnética: -μ * dBz/dz
        Fm = -mag * dBz_dz;

        % Fuerza de fricción: proporcional a la velocidad
        Ff = -gamma * v;

        % Fuerza total: magnética + fricción + peso
        F = Fm + Ff - m * 9.81;

        % Aceleración resultante: F = m*a
        a = F / m;
    end

    % FUNCIÓN ANIDADA: Paso de integración con Runge-Kutta 4
    function [z_next, v_next] = rk4_step(z, v, dt, a_func)
        % rk4_step: realiza un paso de Runge-Kutta 4 para actualizar posición y velocidad.
        % a_func: función anónima que devuelve aceleración en función de (z, v)

        % Primer cálculo (k1)
        k1z = v;
        k1v = a_func(z, v);

        % Segundo cálculo (k2)
        k2z = v + 0.5 * dt * k1v;
        k2v = a_func(z + 0.5 * dt * k1z, v + 0.5 * dt * k1v);

        % Tercer cálculo (k3)
        k3z = v + 0.5 * dt * k2v;
        k3v = a_func(z + 0.5 * dt * k2z, v + 0.5 * dt * k2v);

        % Cuarto cálculo (k4)
        k4z = v + dt * k3v;
        k4v = a_func(z + dt * k3z, v + dt * k3v);

        % Cálculo final de z y v promediando las pendientes
        z_next = z + (dt / 6) * (k1z + 2*k2z + 2*k3z + k4z);
        v_next = v + (dt / 6) * (k1v + 2*k2v + 2*k3v + k4v);
    end

end
