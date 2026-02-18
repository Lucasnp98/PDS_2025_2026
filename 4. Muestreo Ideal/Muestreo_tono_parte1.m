%% MUESTREO IDEAL DE UN COSENO - VERSION LIMPIA (Hz o frecuencia digital)
% 1) Tiempo: coseno continuo + coseno muestreado
% 2) Espectro baseband del coseno muestreado
% 3) Réplicas del espectro
%
% TODO en magnitud lineal
% TODO en el MISMO AZUL (sin ciclo de colores)
% SIN textos ni marcadores

clear; close all; clc;

%% =======================
% Parámetros
%% =======================
fs   = 200;        % Hz
Ts   = 1/fs;
A    = 1;          % amplitud
f0   = 50;         % Hz
Krep = 1;          % nº de réplicas a cada lado

Tobs = 2.0;        % s
Nfft = 8192;       % puntos FFT

% -------- SELECTOR DE EJE DE FRECUENCIA --------
use_digital_freq = true;     % false -> Hz | true -> frecuencia digital
digital_mode     = "w";     % "nu" -> ν=f/fs | "w" -> ω (rad/muestra)
% ----------------------------------------------

% Color único (azul)
blue = [0 0.4470 0.7410];

%% =========================================================
% 1) TIEMPO: coseno continuo vs muestreado
%% =========================================================
dt = 1e-4;
t  = 0:dt:Tobs;
x  = A*cos(2*pi*f0*t);

n   = 0:floor(Tobs*fs);
t_s = n*Ts;
x_n = A*cos(2*pi*f0*t_s);

figure('Name','(1) Tiempo'); grid on; hold on;
plot(t, x, 'LineWidth', 1.8, 'Color', blue);
stem(t_s, x_n, 'filled', 'Color', blue, 'MarkerFaceColor', blue);
xlim([0 10*Ts]);
xlabel('Tiempo (s)');
ylabel('Amplitud');
title(sprintf('Coseno continuo vs muestreado | fs=%.2f Hz, f0=%.2f Hz', fs, f0));
legend('x(t)', 'x[n]', 'Location','best');

%% =========================================================
% 2) ESPECTRO BASEBAND
%% =========================================================
X = fftshift(fft(x_n, Nfft));
f = (-Nfft/2:Nfft/2-1) * (fs/Nfft);
mag = abs(X);

if ~use_digital_freq
    fx = f;
    xlab = 'Frecuencia (Hz)';
    xlim_bb = [-fs/2 fs/2];
    nyq = fs/2;
else
    if digital_mode == "nu"
        fx = f/fs;
        xlab = '\nu = f/f_s (ciclos/muestra)';
        xlim_bb = [-0.5 0.5];
        nyq = 0.5;
    else
        fx = 2*pi*(f/fs);
        xlab = '\omega (rad/muestra)';
        xlim_bb = [-pi pi];
        nyq = pi;
    end
end

figure('Name','(2) Espectro baseband'); grid on; hold on;
plot(fx, mag, 'LineWidth', 1.8, 'Color', blue);
xlim(xlim_bb);
xlabel(xlab);
ylabel('Magnitud');
title('Espectro del coseno muestreado (banda base)');

% Si quieres que estas líneas también sean azules, déjalas así:
xline(-nyq, '--', 'Color', blue);
xline( nyq, '--', 'Color', blue);

%% =========================================================
% 3) RÉPLICAS DEL ESPECTRO (TODAS EN AZUL)
%% =========================================================
figure('Name','(3) Réplicas'); grid on; hold on;

for k = -Krep:Krep
    if ~use_digital_freq
        fx_rep = f + k*fs;
    else
        if digital_mode == "nu"
            fx_rep = (f/fs) + k;
        else
            fx_rep = 2*pi*(f/fs + k);
        end
    end
    plot(fx_rep, mag, 'LineWidth', 1.5, 'Color', blue);
end

if ~use_digital_freq
    xlim([-(Krep+0.5)*fs (Krep+0.5)*fs]);
    for m = -Krep:Krep
        xline(m*fs, ':', 'Color', blue);
    end
else
    if digital_mode == "nu"
        xlim([-(Krep+0.5) (Krep+0.5)]);
        for m = -Krep:Krep
            xline(m, ':', 'Color', blue);
        end
    else
        xlim([-(Krep+0.5)*2*pi (Krep+0.5)*2*pi]);
        for m = -Krep:Krep
            xline(2*pi*m, ':', 'Color', blue);
        end
    end
end

xlabel(xlab);
ylabel('Magnitud');
title(sprintf('Réplicas del espectro del coseno (Krep=%d)', Krep));
