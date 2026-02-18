%% Convolucion_discreta_paso_a_paso.m
% Muestra paso a paso la convolución y[n] = x[n] * h[n]
% - Primero muestra h(n)
% - Luego muestra h volteada h(-n)
% - Luego desplaza h(-n) de izquierda a derecha, y va dibujando y[n] "en tiempo real"
% - Control manual: botones "Siguiente" y "Voltear/Continuar" + tecla Espacio/Enter

clear; close all; clc;

%% ====== 1) Define tus secuencias (EDITA AQUI) ======
x  = [1 1 1 1 1];      nx0 = 0;      % x[n] empieza en n = nx0
h  = [1 2 3];         nh0 = 0;      % h[n] empieza en n = nh0

% (Opcional) velocidad si quisieras modo "auto"
auto_mode = false;    % true = avanza solo; false = manual
auto_pause = 0.25;    % segundos entre pasos en auto

%% ====== 2) Construye ejes n para x y h ======
nx = nx0 : nx0 + numel(x) - 1;
nh = nh0 : nh0 + numel(h) - 1;

%% ====== 3) Convolución "real" para referencia (la iremos dibujando) ======
y = conv(x, h);
ny = (nx(1) + nh(1)) : (nx(end) + nh(end));

%% ====== 4) Preparación figura y UI ======
fig = figure('Name','Convolución discreta paso a paso','Color','w', ...
    'NumberTitle','off','KeyPressFcn',@(~,~) setappdata(gcf,'go',true));

% Appdata para control de pasos
setappdata(fig,'go',false);
setappdata(fig,'flip_ok',false);

% Botón: siguiente
uicontrol('Style','pushbutton','String','Siguiente (o ESPACIO/ENTER)', ...
    'Units','normalized','Position',[0.70 0.94 0.28 0.05], ...
    'FontSize',10,'Callback',@(~,~) setappdata(fig,'go',true));

% Botón: voltear / continuar
uicontrol('Style','pushbutton','String','Voltear/Continuar', ...
    'Units','normalized','Position',[0.02 0.94 0.25 0.05], ...
    'FontSize',10,'Callback',@(~,~) setappdata(fig,'flip_ok',true));

% Subplots
ax1 = subplot(3,1,1);  % x[n]
ax2 = subplot(3,1,2);  % h[n] / h[-n] desplazada
ax3 = subplot(3,1,3);  % y[n] construyéndose

% Rango global para que no "salten" los ejes
nmin = min([nx, nh, ny]);
nmax = max([nx, nh, ny]);

% y límites verticales razonables
ymax = max([abs(x(:)); abs(h(:)); abs(y(:)); 1]);
yL = ymax + 0.5;

%% ====== 5) Dibuja x[n] fijo arriba ======
axes(ax1); cla(ax1);

yLx = max(abs(x)) + 0.5;
stem(ax1, nx, x, 'filled','LineWidth',1.5); grid(ax1,'on');
xlim(ax1,[nmin-1, nmax+1]); ylim(ax1,[-yLx yLx]);
title(ax1,'x[n]');
xlabel(ax1,'n'); ylabel(ax1,'Amplitud');

%% ====== 6) PASO 1: muestra h[n] ======
axes(ax2); cla(ax2);
stem(ax2, nh, h, 'filled','LineWidth',1.5); grid(ax2,'on');
yLh = max(abs(h)) + 1;
xlim(ax2,[nmin-1, nmax+1]); ylim(ax2,[-yLh yLh]);
title(ax2,'h[n] (original)');
xlabel(ax2,'n'); ylabel(ax2,'Amplitud');

% y[n] aún vacío
axes(ax3); cla(ax3);
stem(ax3, ny, zeros(size(ny)), 'LineWidth',1.5); grid(ax3,'on');
yLy = max(abs(y)) + 0.5;
xlim(ax3,[nmin-1, nmax+1]); ylim(ax3,[-yLy yLy]);
title(ax3,'y[n] = x[n] * h[n] (se va construyendo)');
xlabel(ax3,'n'); ylabel(ax3,'Amplitud');

drawnow;

% Espera a que el profe quiera voltear
wait_for_flip(fig, auto_mode, auto_pause);

%% ====== 7) PASO 2: voltear h -> h[-n] ======
% h volteada en muestras, y sus índices cambian: si h estaba en nh, h_flip está en -nh
h_flip = fliplr(h);
nh_flip = -fliplr(nh);    % esto equivale a -(nh) pero alineado con h_flip

axes(ax2); cla(ax2);
stem(ax2, nh_flip, h_flip, 'filled','LineWidth',1.5); grid(ax2,'on');
xlim(ax2,[nmin-1, nmax+1]); ylim(ax2,[-yLh yLh]);
title(ax2,'h[-n] (volteada)');
xlabel(ax2,'n'); ylabel(ax2,'Amplitud');
drawnow;

% Espera a que el profe quiera empezar a desplazar
wait_for_next(fig, auto_mode, auto_pause);

%% ====== 8) PASO 3: desplazar h[-n] y construir y[n] ======
y_build = nan(size(ny));  % lo iremos rellenando

% En convolución: y[n] = sum_k x[k] h[n-k]
% Equivalente con "h[-n] desplazada": en el paso para un n0, alineas h_flip en índices (n0 - nx)
% Más simple: recorremos n0 = ny(1) ... ny(end), calculamos suma, y mostramos h_flip desplazada a n0.

for idx = 1:numel(ny)
    n0 = ny(idx);

    % h_flip desplazada: h[n0-k] se ve como h_flip en índices (n0 - nx) cuando lo superpones con x[nx]
    % Para visualizar: los índices de h_flip desplazada serán:
    nh_shift = n0 - nx;       % lugares donde "cae" h[n0-k] cuando k recorre nx
    h_shift = zeros(size(nx));

    % Para cada k=nx(i), necesitamos h[n0-k]. Vamos a buscarlo en (nh, h)
    % Mapeo directo:
    for i = 1:numel(nx)
        kk = nx(i);
        target_n = n0 - kk;   % índice dentro de h[n]
        j = find(nh == target_n, 1);
        if ~isempty(j)
            h_shift(i) = h(j);
        else
            h_shift(i) = 0;
        end
    end

    % Suma de productos
    y_build(idx) = sum(x .* h_shift);

    % --- Plot de h volteada y desplazada ---
    axes(ax2); cla(ax2);
    % Dibuja x en gris claro como referencia visual en el subplot central
    % stem(ax2, nx, x, 'LineWidth',1.0); hold(ax2,'on');
    % Dibuja h desplazada "en los puntos de k" (misma rejilla que x) para ver el solapamiento
    stem(ax2, nx, h_shift, 'filled','LineWidth',1.5);
    yline(ax2,0,'-');
    grid(ax2,'on');
    xlim(ax2,[nmin-1, nmax+1]); ylim(ax2,[-yLh yLh]);
    title(ax2, sprintf('Paso n = %d: h[n-k] superpuesta sobre x[k] (control manual)', n0));
    xlabel(ax2,'n'); ylabel(ax2,'Amplitud');
    legend(ax2, {'x[k]','h[n-k]'}, 'Location','northeast');

    % --- Plot de y[n] construido hasta ahora ---
    axes(ax3); cla(ax3);
    stem(ax3, ny, zeros(size(ny)),'LineWidth',1.0); hold(ax3,'on');
    stem(ax3, ny(1:idx), y_build(1:idx), 'filled','LineWidth',1.5);
    xline(ax3,n0,'--');
    grid(ax3,'on');
    xlim(ax3,[nmin-1, nmax+1]); ylim(ax3,[-yLy yLy]);
    title(ax3,'y[n] construido (hasta el paso actual)');
    xlabel(ax3,'n'); ylabel(ax3,'Amplitud');

    drawnow;

    % Espera control manual (o auto)
    wait_for_next(fig, auto_mode, auto_pause);
end

%% ====== 9) Fin: muestra resultado final y comparación ======
axes(ax3); cla(ax3);
stem(ax3, ny, y, 'filled','LineWidth',1.5); grid(ax3,'on');
xlim(ax3,[nmin-1, nmax+1]); ylim(ax3,[-yL yL]);
title(ax3,'Resultado final: y[n] = conv(x,h)');
xlabel(ax3,'n'); ylabel(ax3,'Amplitud');

disp('Listo. Convolución final dibujada.');

%% ====== FUNCIONES DE CONTROL ======
function wait_for_next(fig, auto_mode, auto_pause)
    if auto_mode
        pause(auto_pause);
        return;
    end
    setappdata(fig,'go',false);
    % Espera a botón o tecla (KeyPressFcn pone go=true)
    while ishandle(fig) && ~getappdata(fig,'go')
        pause(0.01);
    end
end

function wait_for_flip(fig, auto_mode, auto_pause)
    if auto_mode
        pause(auto_pause);
        return;
    end
    setappdata(fig,'flip_ok',false);
    % Espera al botón "Voltear/Continuar"
    while ishandle(fig) && ~getappdata(fig,'flip_ok')
        pause(0.01);
    end
end
