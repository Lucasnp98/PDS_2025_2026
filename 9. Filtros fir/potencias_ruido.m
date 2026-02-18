%% Ruido blanco gaussiano -> FIR notch en fd = ±0.25 -> Potencias
clear; close all; clc;

% Parámetros
N     = 200000;   % nº muestras
sigma = 0.5;      % desviación típica

% 1) Ruido blanco gaussiano
x = sigma * randn(N,1);

% 2) FIR notch en fd = ±0.25  (ceros en z = ±j)
h = [0.5 0 0.5];

% 3) Filtrar
y = filter(h, 1, x);

% 4) Potencias (definición)
P_in_definicion  = mean(x.^2);
P_out_definicion = mean(y.^2);

fprintf('P_in  = %.10f\n', P_in_definicion);
fprintf('P_out = %.10f\n', P_out_definicion);

%% Autocorrelación (biased => Rx(0)=mean(x^2))
[Rx_in, lags] = xcorr(x, 'biased');
[Rx_out, ~]   = xcorr(y, 'biased');

idx0 = find(lags==0, 1, 'first');   % robusto

Rx0 = Rx_in(idx0);
Ry0 = Rx_out(idx0);

fprintf('\nPotencia desde autocorrelación:\n');
fprintf('Rx_in(0)  = %.10f\n', Rx0);
fprintf('Rx_out(0) = %.10f\n', Ry0);






%% ===============================================================
%% Graficos de las funciones de autocorrelacion
%% ================================================================



%% ============================================================
%% PRINT BONITO: COMPARACIÓN DE POTENCIAS
%% ============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('COMPARACIÓN DE POTENCIA DEL RUIDO\n');
fprintf('============================================================\n');

fprintf('\n--- POTENCIA A LA ENTRADA ---\n');
fprintf('Desde definición:        P_in  = mean(x[n]^2)  = %.10f\n', P_in_definicion);
fprintf('Desde autocorrelación:   Rx_in(0)              = %.10f\n', Rx0);
fprintf('Diferencia absoluta:     |error|               = %.3e\n', abs(P_in_definicion - Rx0));

fprintf('\n--- POTENCIA A LA SALIDA ---\n');
fprintf('Desde definición:        P_out = mean(y[n]^2)  = %.10f\n', P_out_definicion);
fprintf('Desde autocorrelación:   Rx_out(0)             = %.10f\n', Ry0);
fprintf('Diferencia absoluta:     |error|               = %.3e\n', abs(P_out_definicion - Ry0));

fprintf('\n============================================================\n');
fprintf('CONCLUSIÓN:\n');
fprintf('La potencia coincide con el valor de la autocorrelación en\n');
fprintf('el origen:  P = Rx(0)\n');
fprintf('============================================================\n\n');
