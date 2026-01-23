%% sinc(f*tau) en Hz + ceros (rojo) + "deltas" en múltiplos de fs
clear; close all; clc;

%% ======================
% PARÁMETROS (CAMBIA SOLO ESTO)
% ======================
Ts  = 1e-3;          % periodo de muestreo (s)
tau = Ts/3;         % ancho del pulso (s)  <-- prueba Ts/3, Ts/4, etc.

%% ======================
% DERIVADOS
% ======================
fs = 1/Ts;           % frecuencia de muestreo (Hz)

%% ======================
% EJE DE FRECUENCIA (Hz)
% ======================
fMax = 8000;                          % Hz (ajusta si cambias Ts)
f = linspace(-fMax, fMax, 4001);

%% ======================
% SINC CONTINUA
% MATLAB: sinc(x)=sin(pi*x)/(pi*x)
% ======================
S = abs( sinc(f * tau) );

%% ======================
% CEROS DE LA SINC: f = k/tau
% ======================
kMax = floor(fMax * tau);
k = 1:kMax;
fZeros = k / tau;

%% ======================
% "DELTAS" EN MÚLTIPLOS DE fs: f = n*fs
% ======================
nMax = floor(fMax / fs);
n = -nMax:nMax;
fDeltas = n * fs;

% Altura visual de los deltas (solo para dibujar)
deltaHeight = 0.35;

%% ======================
% GRÁFICA
% ======================
figure;

% sinc continua
plot(f, S, 'b', 'LineWidth', 1.6); grid on; hold on;

% deltas (peine) en múltiplos de fs
stem(fDeltas, deltaHeight*ones(size(fDeltas)), 'k', 'filled', 'LineWidth', 1.0);

% ceros de la sinc (rojo)
stem([fZeros, -fZeros], 0.05*ones(1, 2*numel(fZeros)), 'r', 'filled', 'LineWidth', 1.0);

xlabel('Frecuencia (Hz)');
ylabel('Magnitud (escala relativa)');
title(sprintf('|sinc(f\\tau)| con deltas en n f_s   (Ts = %.2f ms,  \\tau = %.2f ms)', ...
               Ts*1e3, tau*1e3));

legend('|sinc(f\tau)|','\delta en n f_s','ceros de sinc','Location','best');

%% ======================
% TEXTO CLAVE
% ======================
fprintf('Ts = %.3g s  -> fs = %.0f Hz\n', Ts, fs);
fprintf('tau = %.3g s  -> ceros en f = k/tau = %.0f, %.0f, ... Hz\n', ...
        tau, 1/tau, 2/tau);
fprintf('Deltas dibujadas en f = n*fs dentro de ±%.0f Hz.\n', fMax);
