%% ============================================================
%% SNR a la entrada y a la salida de un FIR (senal + ruido)
%% ============================================================

clear; close all; clc;

% Parametros
N     = 200000;
n     = (0:N-1).';
sigma = 1;                 % std del ruido

% Senal limpia (sin ruido)
A1 = 2;   fd1 = 0.15;         % tono "bueno"
A2 = 1;   fd2 = 0.25;         % tono que el FIR anula

senal_sin_ruido = A1*cos(2*pi*fd1*n) + A2*cos(2*pi*fd2*n);

% Ruido blanco gaussiano
ruido = sigma * randn(N,1);

% Senal observada (senal + ruido)
senal_con_ruido = senal_sin_ruido + ruido;

% FIR notch en fd = ±0.25
h = [0.5 0 0.5];

% Salida del filtro
salida = filter(h,1,senal_con_ruido);

% Filtrar por separado senal y ruido
salida_senal = filter(h,1,senal_sin_ruido);
salida_ruido = filter(h,1,ruido);

%% Potencias
Ps_in  = mean(senal_sin_ruido.^2);
Pw_in  = mean(ruido.^2);

Ps_out = mean(salida_senal.^2);
Pw_out = mean(salida_ruido.^2);

%% SNR
SNR_in  = Ps_in / Pw_in;
SNR_out = Ps_out / Pw_out;

SNR_in_dB  = 10*log10(SNR_in);
SNR_out_dB = 10*log10(SNR_out);

%% Prints
fprintf('\n============================================================\n');
fprintf('SNR A LA ENTRADA Y A LA SALIDA (FIR)\n');
fprintf('============================================================\n');

fprintf('\n--- ENTRADA ---\n');
fprintf('P_senal,in  = %.10f\n', Ps_in);
fprintf('P_ruido,in  = %.10f\n', Pw_in);
fprintf('SNR_in      = %.6f  (%.3f dB)\n', SNR_in, SNR_in_dB);

fprintf('\n--- SALIDA ---\n');
fprintf('P_senal,out = %.10f\n', Ps_out);
fprintf('P_ruido,out = %.10f\n', Pw_out);
fprintf('SNR_out     = %.6f  (%.3f dB)\n', SNR_out, SNR_out_dB);

fprintf('\n--- MEJORA ---\n');
fprintf('Delta SNR   = %.3f dB\n', SNR_out_dB - SNR_in_dB);

fprintf('============================================================\n\n');





%% ============================================================
%% PLOTS EN EL TIEMPO Y EN FRECUENCIA (senal + ruido)
%% ============================================================

% -------- PARAMETROS FFT --------
Nfft = 8192;
fd = (-Nfft/2:Nfft/2-1)/Nfft;

% FFT entrada y salida
X = fftshift(fft(senal_con_ruido, Nfft))/N;
Y = fftshift(fft(salida, Nfft))/N;

Xmag = abs(X);
Ymag = abs(Y);

%% ============================================================
%% PLOT EN EL TIEMPO
%% ============================================================

Ns = 500; % numero de muestras a mostrar

figure;

subplot(2,1,1)
plot(n(1:Ns), senal_con_ruido(1:Ns), 'LineWidth', 1.5);
grid on;
title('Senal + ruido a la entrada');
xlabel('n');
ylabel('Amplitud');

subplot(2,1,2)
plot(n(1:Ns), salida(1:Ns), 'LineWidth', 1.5);
grid on;
title('Senal + ruido a la salida del FIR');
xlabel('n');
ylabel('Amplitud');


%% ============================================================
%% PLOT EN FRECUENCIA
%% ============================================================

figure;

subplot(2,1,1)
plot(fd, Xmag, 'LineWidth', 1.5);
grid on;
title('Espectro a la entrada');
xlabel('f_d (ciclos/muestra)');
ylabel('|X(f_d)|');
xlim([-0.5 0.5]);

subplot(2,1,2)
plot(fd, Ymag, 'LineWidth', 1.5);
grid on;
title('Espectro a la salida del FIR');
xlabel('f_d (ciclos/muestra)');
ylabel('|Y(f_d)|');
xlim([-0.5 0.5]);


x = 5 + randn(1000000,1);

P = mean(x.^2);

[Rx,lags] = xcorr(x,'biased');
idx0 = find(lags==0,1);
Rx0 = Rx(idx0);

fprintf("P = %.4f\n",P)
fprintf("Rx(0) = %.4f\n",Rx0)
