%% Pulso -> repetir con FOR -> pinchos en frecuencia
numero_periodos = 10;     % 1, 2, 5, 20, 100...

clearvars -except numero_periodos
close all; clc;

%% Parámetros sencillos
A   = 1;
tau = 0.2e-3;     % duración del pulso
T0  = 1e-3;       % periodo
Fs  = 2e6;        % muestreo

dt = 1/Fs;
f0 = 1/T0;

%% Número de muestras (sin floor ni pollas)
muestras_pulso   = tau * Fs;
muestras_periodo = T0  * Fs;

%% 1️⃣ Definir UN pulso (unos y ceros)
pulso = [ A*ones(1, muestras_pulso), zeros(1, muestras_periodo - muestras_pulso) ];

%% 2️⃣ Repetir el pulso CON UN FOR
x = [];

for k = 1:numero_periodos
    x = [x pulso];   % 👈 concatenar, lento pero CLARO
end

%% 3️⃣ Eje temporal centrado en 0
N = length(x);
t = (0:N-1)*dt;
t = t - mean(t);   % centrar sin pensar

%% 4️⃣ FFT
X = fftshift(fft(x));
f = (-N/2:N/2-1)*(Fs/N);
X = dt * X;

%% 5️⃣ Dibujos
figure('Color','w');

subplot(2,1,1)
plot(t*1e3, x, 'LineWidth', 1.5);
grid on;
xlabel('t (ms)')
ylabel('x(t)')
title(['Pulso repetido ', num2str(numero_periodos), ' veces (centrado)'])

subplot(2,1,2)
plot(f/1e3, abs(X), 'LineWidth', 1.5);
grid on;
xlabel('f (kHz)')
ylabel('|X(f)|')
title('Espectro: aparecen pinchos al aumentar N')
xlim([-200 200])
