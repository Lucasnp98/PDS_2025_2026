%% demo_notch_general_dos_plots.m
clear; close all; clc;

%% ==============================
%% 1) SEÑAL GENERAL
%% ==============================

N = 1000;
n = (0:N-1)';

fd_list = [0.15 0.3 ];    % modificable
A_list  = [2  1 ];       % modificable

if numel(fd_list) ~= numel(A_list)
    error('fd_list y A_list deben tener la misma longitud');
end

x = zeros(size(n));

for k = 1:numel(fd_list)
    x = x + A_list(k)*cos(2*pi*fd_list(k)*n);
end


%% ==============================
%% 2) FILTRO FIR (modificable)
%% ==============================

h = [1/2 0 1/2];


%% ==============================
%% 3) SALIDA Y POTENCIAS
%% ==============================

y = filter(h,1,x);

Px = mean(x.^2);
Py = mean(y.^2);

fprintf('\nPotencias:\n');
fprintf('Entrada Px = %.6f\n', Px);
fprintf('Salida  Py = %.6f\n', Py);


%% ==============================
%% 4) FFTs
%% ==============================

Nfft = 4096;
fd = (-Nfft/2:Nfft/2-1)/Nfft;

X = fftshift(fft(x,Nfft))/N;
Y = fftshift(fft(y,Nfft))/N;

Xmag = abs(X);
Ymag = abs(Y);

w = 2*pi*fd;
H = freqz(h,1,w);
Hmag = abs(H);


%% ==============================
%% 5) FRECUENCIAS AUTOMÁTICAS
%% ==============================

ceros = roots(h);
fd_ceros = unique(round(angle(ceros)/(2*pi),6));











%% CÁLCULO DE LAS RESPUESTAS EN FRECUENCIA:



%% ==============================
%% EXTRA) H(e^{jw}) en las frecuencias del/los cosenos (ROBUSTO)
%% ==============================
L = numel(h);
m = (0:L-1).';                 % índices del FIR
wd_list = 2*pi*fd_list(:);     % rad/muestra

% Hk(k) = sum_m h[m]*exp(-j*wd_k*m)
Hk = (h(:).' * exp(-1j*(m*wd_list.'))).';   % tamaño: [K x 1]

Ak_in = A_list(:);
Ak_out_est = Ak_in .* abs(Hk);             % amplitud estimada salida

fprintf('\n--- H en las frecuencias de la señal ---\n');
for k = 1:numel(fd_list)
    fprintf('k=%d, fd=%.6f -> |H|=%.6f, ∠H=%.2f deg | A_in=%.3f -> A_out≈%.3f\n', ...
        k, fd_list(k), abs(Hk(k)), angle(Hk(k))*180/pi, Ak_in(k), Ak_out_est(k));
end

%% ==============================
%% EXTRA) H en las frecuencias de la(s) componente(s) del coseno
%% ==============================




%% ==============================
%% 6) PLOT 1: ENTRADA vs FILTRO
%% ==============================

figure;

yyaxis left
plot(fd, Xmag, 'LineWidth',1.5)
ylabel('|X(f_d)|')

yyaxis right
plot(fd, Hmag, 'LineWidth',1.5)
ylabel('|H(f_d)|')

xlabel('f_d')
title('Espectro de entrada y respuesta del filtro')
grid on
xlim([-0.5 0.5])
hold on

% marcar frecuencias señal
for k = 1:numel(fd_list)
    xline(fd_list(k),'--');
    xline(-fd_list(k),'--');
end

% marcar notches filtro
for k = 1:numel(fd_ceros)
    xline(fd_ceros(k),'--r','LineWidth',1.5);
end

legend('Entrada |X|','Filtro |H|')












%% ==============================
%% 7) PLOT 2: SOLO ESPECTRO SALIDA
%% ==============================

figure;

plot(fd, Ymag, 'LineWidth',1.5)
xlabel('f_d')
ylabel('|Y(f_d)|')

title('Espectro de la salida del filtro')
grid on
xlim([-0.5 0.5])
hold on

% marcar notches automáticamente
for k = 1:numel(fd_ceros)
    xline(fd_ceros(k),'--r','LineWidth',1.5);
end
