function animar_con_fondo(zm, imagen, nombre_video)
% Esta función genera una animación de la trayectoria de un objeto (por ejemplo, un imán)
% cayendo a lo largo del eje z, con una imagen de fondo representando el campo magnético.
%
% Parámetros:
% zm           -> Vector de posiciones en z del objeto a animar
% imagen       -> Archivo de imagen (PNG) que se usa como fondo
% nombre_video -> Nombre del archivo de salida en video (MP4 o AVI)

    % Define el rango físico de la imagen en el plano X-Z
    % Esto debe coincidir con el rango que usaste al graficar el campo
    x_range = [-6, 6];  % Rango en el eje X
    z_range = [-6, 6];  % Rango en el eje Z

    % Crear el objeto para escribir video
    % Puedes usar 'MPEG-4' o 'Motion JPEG AVI' (más compatible)
    writerObj = VideoWriter(nombre_video, 'Motion JPEG AVI');

    % Configura la velocidad de la animación (fotogramas por segundo)
    writerObj.FrameRate = 30;  % Puedes poner un valor menor para animaciones más lentas
    open(writerObj);           % Abre el archivo de video para escritura

    % Cargar imagen de fondo (campo magnético u otra ilustración)
    img = imread(imagen);

    % Crear ventana de figura para animación
    figure(100); clf  % Nueva figura con ID 100, y limpiar contenido previo

    % Dibujar la imagen de fondo invertida verticalmente (flipud)
    imagesc(x_range, z_range, flipud(img));
    axis xy;          % Asegura que el eje z crezca hacia arriba (sistema cartesiano)
    hold on;
    xlabel('x (m)');
    ylabel('z (m)');
    title('Trayectoria sobre campo magnético');

    % Animar punto por punto según las posiciones en zm
    for i = 1:length(zm)
        cla;  % Borra el contenido anterior sin cerrar la figura

        % Redibujar fondo en cada frame (porque cla borra todo)
        imagesc(x_range, z_range, flipud(img));
        axis xy;
        hold on;

        % Dibuja el objeto en su posición actual (0 en x, zm(i) en z)
        scatter(0, zm(i), 100, 'r', 'filled');  % Punto rojo grande

        % Limita los ejes para que no cambien durante la animación
        xlim(x_range);
        ylim(z_range);

        % Captura la figura actual como un frame
        frame = getframe(gcf);

        % Escribe el frame en el archivo de video
        writeVideo(writerObj, frame);
    end

    % Finaliza el archivo de video
    close(writerObj);

    % Mensaje de confirmación en consola
    disp(['Video guardado como ', nombre_video]);
end
