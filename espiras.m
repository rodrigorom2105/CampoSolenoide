function [Px, Py, Pz, dx, dy, dz] = espiras(nl, N, R, sz)
% Esta función genera y visualiza un conjunto de 'nl' espiras circulares
% con 'N' puntos por espira, de radio 'R' y separadas entre sí por 'sz' en el eje z.
% También calcula la dirección de la corriente en cada segmento de espira.

% -------------------------
% PARTE 1: GENERACIÓN DE ESPIRAS
% -------------------------

    % Paso angular entre puntos (división uniforme del círculo)
    dtheta = 2*pi / N;

    % Vector de ángulos para una espira completa
    ang = 0:dtheta:(2*pi - dtheta);  % Se generan N puntos desde 0 hasta 2π (sin repetir el 0 final)

    % Inicializamos la variable de índice para colocación de datos
    s = 1;

    % Bucle que itera por cada espira (nl en total)
    for i = 1:nl
        % Coordenadas X de los puntos de la espira actual
        Px(s:s+N-1) = R * cos(ang);

        % Coordenadas Y de los puntos de la espira actual
        Py(s:s+N-1) = R * sin(ang);

        % Coordenadas Z (alturas) para los puntos de esta espira
        % Cada espira se coloca a una altura diferente, comenzando desde -nl/2*sz hasta +nl/2*sz
        Pz(s:s+N-1) = -nl/2*sz + (i-1)*sz;

        % Componentes del vector de dirección de corriente (tangente al círculo en cada punto)
        % La corriente en una espira circular tiene dirección tangente a la trayectoria:
        dx(s:s+N-1) = -Py(s:s+N-1) * dtheta;  % -sin(θ)
        dy(s:s+N-1) =  Px(s:s+N-1) * dtheta;  % cos(θ)

        % Se actualiza el índice para los siguientes N puntos
        s = s + N;
    end

    % Como la corriente fluye en el plano XY, no hay componente en Z
    dz = zeros(1, N*nl);

% -------------------------
% PARTE 2: VISUALIZACIÓN
% -------------------------

    % Crea una nueva figura para visualizar las espiras
    figure(1);

    % Dibuja vectores en 3D representando la dirección de la corriente
    % en cada punto (posición: Px,Py,Pz - dirección: dx,dy,dz)
    quiver3(Px, Py, Pz, dx, dy, dz, 0.5, '-r', 'LineWidth', 2);

    % Ajusta el ángulo de visualización 3D
    view(-34, 33);

    % Título y etiquetas de los ejes
    title('Corriente en espiras');
    xlabel('x'); ylabel('y'); zlabel('z');

    % Fuerza a que los ejes tengan la misma escala (importante para formas circulares)
    axis equal;
end