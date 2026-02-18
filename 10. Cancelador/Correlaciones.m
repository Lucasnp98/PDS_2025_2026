%% ============================================================
%% Autocorrelación entrada/salida y potencia = Rx(0)
%% ============================================================

clear; close all; clc;

% Parámetros
N     = 1000000;
sigma = 1;

% Ruido blanco gaussiano
x = sigma * randn(N,1);

% FIR notch en fd = ±0.25
h = [0.5 0 0.5];
y = filter(h, 1, x);

%% 1) Potencia por definición
P_in  = mean(x.^2);
P_out = mean(y.^2);

%% 2) Autocorrelaciones (biased => Rx(0) = mean(x^2))
[Rx_in, lags] = xcorr(x, 'biased');
[Rx_out, ~]   = xcorr(y, 'biased');

% Índice del origen (lag = 0)
idx0 = find(lags==0, 1, 'first');

% Valor en el origen (potencia vía autocorrelación)
Rx0_in  = Rx_in(idx0);
Rx0_out = Rx_out(idx0);

%% 3) Prints entrada/salida (definición vs autocorrelación)
fprintf('\n================= POTENCIA = Rx(0) =================\n');
fprintf('Entrada: mean(x^2) = %.10f | Rx_in(0)  = %.10f | error = %.3e\n', ...
        P_in,  Rx0_in,  abs(P_in  - Rx0_in));
fprintf('Salida : mean(y^2) = %.10f | Rx_out(0) = %.10f | error = %.3e\n', ...
        P_out, Rx0_out, abs(P_out - Rx0_out));
fprintf('====================================================\n\n');

%% 4) Plots autocorrelación (zoom alrededor de 0)
L = min([20, idx0-1, length(lags)-idx0]); % evita errores si N es pequeño
rg = (idx0-L):(idx0+L);

figure;
subplot(2,1,1)
stem(lags(rg), Rx_in(rg), 'filled'); grid on;
title('Autocorrelación a la entrada R_{x,in}[k] (zoom)');
xlabel('Lag k'); ylabel('R_{x,in}[k]');

subplot(2,1,2)
stem(lags(rg), Rx_out(rg), 'filled'); grid on;
title('Autocorrelación a la salida R_{x,out}[k] (zoom)');
xlabel('Lag k'); ylabel('R_{x,out}[k]');

% Marcar el valor en el origen
subplot(2,1,1); hold on; stem(0, Rx0_in, 'filled', 'LineWidth', 2);
subplot(2,1,2); hold on; stem(0, Rx0_out,'filled', 'LineWidth', 2);

%% ============================================================
%% COMPARACIÓN DE LAS 3 FORMAS DE CALCULAR LA POTENCIA DE SALIDA
%% ============================================================

% 1) Definición directa
P_out_def = P_out;

% 2) Autocorrelación
P_out_corr = Rx0_out;

% 3) Fórmula teórica del FIR
Pcalculada = P_in * sum(h.^2);

fprintf('\n============================================================\n');
fprintf('POTENCIA DEL RUIDO A LA SALIDA DEL FILTRO FIR\n');
fprintf('============================================================\n');

fprintf('\n1) Por definición directa:\n');
fprintf('   P_out = mean(y[n]^2)        = %.10f\n', P_out_def);

fprintf('\n2) Desde la autocorrelación:\n');
fprintf('   P_out = Rx_out(0)           = %.10f\n', P_out_corr);

fprintf('\n3) Desde la fórmula teórica FIR:\n');
fprintf('   P_out = P_in * sum(h[n]^2)  = %.10f\n', Pcalculada);

fprintf('\n------------------------------------------------------------\n');
fprintf('Errores absolutos:\n');
fprintf('|def - corr| = %.3e\n', abs(P_out_def - P_out_corr));
fprintf('|def - teo | = %.3e\n', abs(P_out_def - Pcalculada));

fprintf('\nCONCLUSIÓN:\n');
fprintf('Las tres formas coinciden: la potencia es única y consistente.\n');
fprintf('============================================================\n\n');






%% FFT de Rx y plot (PSD ~ FFT{Rx})

Nfft = 4096;  

% Reordenar la autocorrelación para que k=0 esté en la posición 1
Rx_shifted = ifftshift(Rx_in);

% FFT de la autocorrelación (Wiener–Khinchin)
Sx = fftshift(fft(Rx_shifted, Nfft));

% Parte real (la imaginaria es error numérico)
Sx = real(Sx);

% Eje de frecuencia normalizada
fd = (-Nfft/2:Nfft/2-1)/Nfft;

figure;
plot(fd, Sx, 'LineWidth', 1.5); grid on;
xlabel('f_d (ciclos/muestra)');
ylabel('S_x(f_d)');
title('Densidad espectral de potencia S_x(f_d) = FFT\{R_x[k]\}');
xlim([-0.5 0.5]);
