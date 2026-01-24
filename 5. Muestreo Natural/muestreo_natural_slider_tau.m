function muestreo_natural_slider_tau()
%% MUESTREO NATURAL (gating) DE UN COSENO - SLIDER PARA TAU (ANCHO DE PULSO)
% (3) incluye superposición de |sinc(f*tau)| (implementada a mano: NO TOOLBOX)

clear; close all; clc;

%% =======================
% Parámetros
%% =======================
fs   = 200;         % Hz
Ts   = 1/fs;

A    = 1;           % amplitud
f0   = 25;          % Hz (tono)

Krep = 3;           % réplicas a cada lado (solo para límites/guías)
Tobs = 10;         % s

Ns_per_Ts = 400;    % puntos por Ts
dt = Ts / Ns_per_Ts;

Nfft = 2^nextpow2(round(Tobs/dt));  % FFT para CTFT aprox

% -------- SELECTOR DE EJE DE FRECUENCIA --------
use_digital_freq = false;     % false -> Hz | true -> frecuencia digital
digital_mode     = 'nu';     % 'nu' -> ν=f/fs | 'w' -> ω (rad/muestra)
% ----------------------------------------------

% Slider tau
tau_min  = 0.01*Ts;
tau_max  = 0.99*Ts;
tau_init = 0.50*Ts;

%% =======================
% Colores (distintos)
%% =======================
col_time_x = [0 0.4470 0.7410];        % azul (x(t))
col_time_y = [0.8500 0.3250 0.0980];   % naranja (y(t))
col_spec_y = [0 0.4470 0.7410];        % azul (|Y|)
col_sinc   = [0.6350 0.0780 0.1840];   % rojo (|sinc|)
col_guides = [0.2 0.2 0.2];            % gris (guías)

%% =======================
% Precomputo ejes
%% =======================
t = 0:dt:Tobs-dt;
x = A*cos(2*pi*f0*t);

f = (-Nfft/2:Nfft/2-1) * (1/(dt*Nfft));  % Hz
[fx_base, xlab, xlim_bb, nyq] = freq_axis(f, fs, use_digital_freq, digital_mode);

%% =======================
% Crear figuras y objetos
%% =======================
% (1) Tiempo
fig1 = figure('Name','(1) Tiempo + Slider \tau (muestreo natural)');
ax1 = axes(fig1); grid(ax1,'on'); hold(ax1,'on');

plot(ax1, t, x, 'LineWidth', 1.8, 'Color', col_time_x);
h_y = plot(ax1, t, zeros(size(t)), '--', 'LineWidth', 1.8, 'Color', col_time_y);

xlabel(ax1, 'Tiempo (s)'); ylabel(ax1, 'Amplitud');
h_title1 = title(ax1, '');
xlim(ax1, [0, 10*Ts]);
legend(ax1, {'x(t)', 'y(t)=x(t)\cdot p(t)'}, 'Location','best');

% (2) Baseband
fig2 = figure('Name','(2) Espectro baseband (muestreo natural)');
ax2 = axes(fig2); grid(ax2,'on'); hold(ax2,'on');

h_base = plot(ax2, fx_base, zeros(size(fx_base)), 'LineWidth', 1.8, 'Color', col_spec_y);
xlim(ax2, xlim_bb);
xlabel(ax2, xlab); ylabel(ax2, 'Magnitud');
title(ax2, 'Espectro baseband de y(t) (muestreo natural)');
xline(ax2, -nyq, '--', 'Color', col_guides);
xline(ax2,  nyq, '--', 'Color', col_guides);

% (3) Amplio (réplicas) + sinc
fig3 = figure('Name','(3) Espectro amplio (réplicas) + sinc - muestreo natural');
ax3 = axes(fig3); grid(ax3,'on'); hold(ax3,'on');

h_wide = plot(ax3, fx_base, zeros(size(fx_base)), 'LineWidth', 1.5, 'Color', col_spec_y);
h_sinc = plot(ax3, fx_base, zeros(size(fx_base)), '--', 'LineWidth', 1.5, 'Color', col_sinc);

legend(ax3, 'off');   
xlabel(ax3, xlab); ylabel(ax3, 'Magnitud');
h_title3 = title(ax3, '');
%legend(ax3, {'|Y(f)| (réplicas)', '|sinc(f\tau)| (envolvente)'}, 'Location','best');

% Límites y guías de réplica
if ~use_digital_freq
    xlim(ax3, [-(Krep+0.5)*fs (Krep+0.5)*fs]);
    for m = -Krep:Krep
        xline(ax3, m*fs, ':', 'Color', col_guides);
    end
else
    if strcmpi(digital_mode,'nu')
        xlim(ax3, [-(Krep+0.5) (Krep+0.5)]);
        for m = -Krep:Krep
            xline(ax3, m, ':', 'Color', col_guides);
        end
    else
        xlim(ax3, [-(Krep+0.5)*2*pi (Krep+0.5)*2*pi]);
        for m = -Krep:Krep
            xline(ax3, 2*pi*m, ':', 'Color', col_guides);
        end
    end
end

%% =======================
% Slider
%% =======================
uicontrol(fig1, 'Style', 'text', 'String', '\tau', ...
    'Units','normalized', 'Position',[0.12 0.02 0.05 0.04], ...
    'BackgroundColor', get(fig1,'Color'), 'FontSize', 11);

h_slider = uicontrol(fig1, 'Style', 'slider', ...
    'Units','normalized', 'Position',[0.18 0.02 0.62 0.04], ...
    'Min', tau_min, 'Max', tau_max, 'Value', tau_init);

h_tau_text = uicontrol(fig1, 'Style', 'text', ...
    'Units','normalized', 'Position',[0.81 0.02 0.17 0.04], ...
    'String', '', 'BackgroundColor', get(fig1,'Color'), 'FontSize', 11);

h_slider.Callback = @(src,~) update_all(get(src,'Value'));

update_all(tau_init);

%% =======================
% Update
%% =======================
    function update_all(tau)
        D = tau / Ts;

        p = mod(t, Ts) < tau;     % gating
        y = x .* p;

        set(h_y, 'YData', y);
        set(h_title1, 'String', sprintf('Muestreo natural | fs=%.1f Hz, f0=%.1f Hz, \\tau=%.4g s (D=%.2f)', fs, f0, tau, D));
        set(h_tau_text, 'String', sprintf('\\tau=%.4g s  (D=%.2f)', tau, D));

        % Espectro (CTFT aprox)
        Y = fftshift(fft(y, Nfft)) * dt;
        MagY = abs(Y);

        set(h_base, 'YData', MagY);
        set(h_wide, 'XData', fx_base, 'YData', MagY);
        set(h_title3, 'String', sprintf('Espectro amplio (réplicas) + sinc | Krep=%d | \\tau=%.4g s (D=%.2f)', Krep, tau, D));

        % sinc propia: sinc(x)=sin(pi*x)/(pi*x)
        S = abs(mysinc(f*tau));
        S = S / (max(S) + eps);   % normaliza para superponer
        set(h_sinc, 'XData', fx_base, 'YData', S);

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
    if strcmpi(digital_mode,'nu')
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
end

%% ===== Helper: sinc propia (sin toolbox) =====
function y = mysinc(x)
y = ones(size(x));
nz = (x ~= 0);
y(nz) = sin(pi*x(nz)) ./ (pi*x(nz));
end
