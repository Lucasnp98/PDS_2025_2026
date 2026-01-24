function muestreo_coseno_slider()
%% MUESTREO IDEAL DE UN COSENO - VERSION CON SLIDER (f0)
% 1) Tiempo: coseno continuo + coseno muestreado
% 2) Espectro baseband del coseno muestreado
% 3) Réplicas del espectro
%
% Slider para variar f0 en tiempo real (actualiza los 3 gráficos).
% TODO en magnitud lineal
% TODO en el MISMO AZUL

clear; close all; clc;

%% =======================
% Parámetros fijos
%% =======================
fs   = 200;        % Hz
Ts   = 1/fs;
A    = 1;          % amplitud
Krep = 1;          % nº de réplicas a cada lado

Tobs = 2.0;        % s
dt   = 1e-4;       % s (para señal continua)
Nfft = 8192;       % puntos FFT

% -------- SELECTOR DE EJE DE FRECUENCIA --------
use_digital_freq = true;     % false -> Hz | true -> frecuencia digital
digital_mode     = "nu";     % "nu" -> ν=f/fs | "w" -> ω (rad/muestra)
% ----------------------------------------------

% Rango del slider (para ver aliasing bien)
f0_min = 0;
f0_max = fs;       % puedes poner fs/2 o fs si quieres ver plegado fuerte
f0_init = 25;

% Color único
blue = [0 0.4470 0.7410];

%% =======================
% Precomputo ejes comunes
%% =======================
t = 0:dt:Tobs;                  % continuo (para plot tiempo)
n = 0:floor(Tobs*fs);           % muestras
t_s = n*Ts;

f = (-Nfft/2:Nfft/2-1) * (fs/Nfft);  % Hz (baseband)
% Conversión eje de frecuencia (según modo)
[fx_base, xlab, xlim_bb, nyq] = freq_axis(f, fs, use_digital_freq, digital_mode);

%% =======================
% Crear figuras y objetos gráficos (actualizables)
%% =======================
% --- (1) Tiempo ---
fig1 = figure('Name','(1) Tiempo + Slider f0');
ax1 = axes(fig1); grid(ax1,'on'); hold(ax1,'on');

h_xt = plot(ax1, t, zeros(size(t)), 'LineWidth', 1.8, 'Color', blue);
h_xn = stem(ax1, t_s, zeros(size(t_s)), 'filled', ...
    'Color', "red", 'MarkerFaceColor','red');

xlim(ax1, [0 10*Ts]);
xlabel(ax1, 'Tiempo (s)'); ylabel(ax1, 'Amplitud');
h_title1 = title(ax1, '');

% --- (2) Espectro baseband ---
fig2 = figure('Name','(2) Espectro baseband');
ax2 = axes(fig2); grid(ax2,'on'); hold(ax2,'on');

h_spec = plot(ax2, fx_base, zeros(size(fx_base)), 'LineWidth', 1.8, 'Color', blue);
xlim(ax2, xlim_bb);
xlabel(ax2, xlab); ylabel(ax2, 'Magnitud');
title(ax2, 'Espectro (banda base)');

% Líneas Nyquist (también azules)
h_nyq1 = xline(ax2, -nyq, '--', 'Color', blue);
h_nyq2 = xline(ax2,  nyq, '--', 'Color', blue);

% --- (3) Réplicas ---
fig3 = figure('Name','(3) Réplicas');
ax3 = axes(fig3); grid(ax3,'on'); hold(ax3,'on');
xlabel(ax3, xlab); ylabel(ax3, 'Magnitud');
h_title3 = title(ax3, '');

% Creamos una línea por réplica (k=-Krep:Krep), todas azules
k_list = -Krep:Krep;
h_rep = gobjects(size(k_list));
for i = 1:numel(k_list)
    h_rep(i) = plot(ax3, nan, nan, 'LineWidth', 1.5, 'Color', blue);
end

% Líneas guía en centros de réplicas (también azules)
h_guides = gobjects(size(k_list));
for i = 1:numel(k_list)
    if ~use_digital_freq
        xc = k_list(i)*fs;
    else
        if digital_mode == "nu"
            xc = k_list(i);
        else
            xc = 2*pi*k_list(i);
        end
    end
    h_guides(i) = xline(ax3, xc, ':', 'Color', blue);
end

% Límites del plot de réplicas
if ~use_digital_freq
    xlim(ax3, [-(Krep+0.5)*fs (Krep+0.5)*fs]);
else
    if digital_mode == "nu"
        xlim(ax3, [-(Krep+0.5) (Krep+0.5)]);
    else
        xlim(ax3, [-(Krep+0.5)*2*pi (Krep+0.5)*2*pi]);
    end
end

%% =======================
% Slider (en la figura 1)
%% =======================
% UI en figura 1 (tiempo) para que sea cómodo
uicontrol(fig1, 'Style', 'text', 'String', 'f0', ...
    'Units','normalized', 'Position',[0.15 0.02 0.05 0.04], ...
    'BackgroundColor', get(fig1,'Color'), 'FontSize', 11);

h_slider = uicontrol(fig1, 'Style', 'slider', ...
    'Units','normalized', 'Position',[0.22 0.02 0.60 0.04], ...
    'Min', f0_min, 'Max', f0_max, 'Value', f0_init);

h_f0_text = uicontrol(fig1, 'Style', 'text', ...
    'Units','normalized', 'Position',[0.83 0.02 0.12 0.04], ...
    'String', sprintf('%.2f Hz', f0_init), ...
    'BackgroundColor', get(fig1,'Color'), 'FontSize', 11);

% Callback: actualizar todo
h_slider.Callback = @(src,~) update_all(src.Value);

% Primera actualización
update_all(f0_init);

%% =======================
% Función de actualización
%% =======================
    function update_all(f0)
        % --- Tiempo ---
        xt = A*cos(2*pi*f0*t);
        xn = A*cos(2*pi*f0*t_s);

        set(h_xt, 'YData', xt);
        set(h_xn, 'YData', xn);

        set(h_title1, 'String', sprintf('Coseno continuo vs muestreado | fs=%.2f Hz, f0=%.2f Hz', fs, f0));
        set(h_f0_text, 'String', sprintf('%.2f Hz', f0));

        % --- Espectro baseband ---
        X = fftshift(fft(xn, Nfft));
        mag = abs(X);

        set(h_spec, 'YData', mag);

        % --- Réplicas: reutilizamos el espectro baseband y lo desplazamos
        for ii = 1:numel(k_list)
            kk = k_list(ii);

            if ~use_digital_freq
                fx_rep = f + kk*fs;
            else
                if digital_mode == "nu"
                    fx_rep = (f/fs) + kk;
                else
                    fx_rep = 2*pi*(f/fs + kk);
                end
            end

            set(h_rep(ii), 'XData', fx_rep, 'YData', mag);
        end

        set(h_title3, 'String', sprintf('Réplicas del espectro (Krep=%d) | f0=%.2f', Krep, f0));
        drawnow limitrate;
    end

end

%% ===== Helper: eje de frecuencia =====
function [fx, xlab, xlim_bb, nyq] = freq_axis(f, fs, use_digital_freq, digital_mode)
if ~use_digital_freq
    fx = f;
    xlab = 'Frecuencia (Hz)';
    xlim_bb = [-fs/2 fs/2];
    nyq = fs/2;
else
    if digital_mode == "nu"
        fx = f/fs;                         % ν
        xlab = '\nu = f/f_s (ciclos/muestra)';
        xlim_bb = [-0.5 0.5];
        nyq = 0.5;
    else
        fx = 2*pi*(f/fs);                  % ω
        xlab = '\omega (rad/muestra)';
        xlim_bb = [-pi pi];
        nyq = pi;
    end
end
end
