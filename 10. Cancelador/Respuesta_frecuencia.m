%% FIR_notch_demo_fullscript.m

% 1) Espectro de entrada (amplitud física) + |H(fd)| en la misma figura
% 2) Espectro de entrada y salida (amplitud física) en la misma figura
% Frecuencia digital fd en [-0.5, 0.5)

clear; close all; clc;

%% Parámetros señal
N = 1000;                     % elige N múltiplo de 10 para que fd=0.3 caiga en bin (0.3*N = entero)
n = (0:N-1)';

fd_list = [0.15 0.3];   % ciclos/muestra
A_list  = [1     1 ];     % amplitudes de cada coseno

x = zeros(size(n));

for k = 1:numel(fd_list)
    x = x + A_list(k)*cos(2*pi*fd_list(k)*n);
end


h = [1, -sqrt(2), 2, -sqrt(2), 1];

%h = [1/2 0 1/2];
% 
%h = [1 2 2 2 1];




% Normaliza a ganancia 1 en DC 
h = h / sum(h);

y = filter(h, 1, x);

%% Rejilla de frecuencia digital completa [-0.5, 0.5)
Nfft = 4096;
fd = (-Nfft/2:Nfft/2-1)/Nfft;    % ciclos/muestra

%% FFTs (amplitud "física")
% Usamos división por N (longitud de la señal en tiempo).
% Para un coseno de amplitud 1 alineado en bin:
%   |X| tendrá picos de 0.5 en +fd y 0.5 en -fd.
X = fftshift(fft(x, Nfft)) / N;
Y = fftshift(fft(y, Nfft)) / N;

Xmag = abs(X);
Ymag = abs(Y);

%% Respuesta en frecuencia del filtro en la misma rejilla
% freqz admite w en rad/muestra (aquí w = 2*pi*fd en [-pi, pi))
w = 2*pi*fd;
H = freqz(h, 1, w);
Hmag = abs(H);

%% FIGURA 1: Espectro de entrada + |H|
figure('Name','Figura 1: Entrada y |H|');
yyaxis left
plot(fd, Xmag, 'LineWidth', 1.5);
ylabel('|X(f_d)| (amplitud por exponencial)');
grid on;

yyaxis right
plot(fd, Hmag, 'LineWidth', 1.5);
ylabel('|H(f_d)| (ganancia)');

xlabel('f_d (ciclos/muestra)');
title('Espectro de entrada (normalizado por N) y respuesta del FIR');
xlim([-0.5 0.5]);


legend('|X|','|H|','Location','best');

%% FIGURA 2: Espectro de entrada y salida
figure('Name','Figura 2: Entrada vs Salida');
plot(fd, Xmag, 'LineWidth', 1.5); hold on;
plot(fd, Ymag, 'LineWidth', 1.5);
grid on;

xlabel('f_d (ciclos/muestra)');
ylabel('Magnitud (normalizada por N)');
title('Espectro de entrada y salida (ambos normalizados por N)');
xlim([-0.5 0.5]);


legend('Entrada |X|','Salida |Y|','Location','best');

%% Extra docente: comprobar ganancia del filtro en las 3 frecuencias
w_test = 2*pi*fd_list;
H_test = freqz(h,1,w_test);

fprintf('\nGanancias del FIR en las frecuencias:\n');
for i = 1:numel(fd_list)
    fprintf('  fd = %.3f  ->  |H| = %.8f\n', fd_list(i), abs(H_test(i)));
end

% Nota:
% - Verás |H| ~ 0 en 0.125 y 0.25 (notches).
% - Verás |H| ~ 1 (aprox) en 0.3 (depende de la normalización por sum(h)).




figure;

%% (1) Plano-z + ceros en rojo
ceros = roots(h);

subplot(1,2,1)
zplane(h,1)
title('Plano-z')
hold on
plot(real(ceros), imag(ceros), 'ro', 'MarkerSize', 10, 'LineWidth', 2)

%% (2) Respuesta en frecuencia + líneas en las frecuencias de los ceros
subplot(1,2,2)
[H,w] = freqz(h,1,4096);
fd = w/(2*pi);

plot(fd, abs(H), 'LineWidth', 1.5)
grid on
title('Respuesta en frecuencia')
xlabel('f_d')
hold on

% --- Frecuencias digitales asociadas a los ceros ---
% ángulos en [-pi, pi)
ang = angle(ceros);
fd_ceros = ang/(2*pi);              % fd en [-0.5, 0.5)

% Nos quedamos con las frecuencias positivas (0..0.5), quitando duplicados
fd_pos = fd_ceros(fd_ceros > 0);
fd_pos = unique(round(fd_pos, 6));  % agrupar numéricamente

% Dibujar líneas verticales en esas frecuencias
for k = 1:numel(fd_pos)
    xline(fd_pos(k), '--r');
end

xlim([0 0.5])
