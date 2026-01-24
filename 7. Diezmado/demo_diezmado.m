%% DEMO DIEZMADO (DECIMATION) CON SINC -> RECT EN FRECUENCIA (SCRIPT COMPLETO)
% Señal continua: x(t) = sinc(2*B*t)
% (MATLAB sinc(u) = sin(pi*u)/(pi*u))
%
% CTFT ideal:
%   X(f) = (1/(2B)) * rect(f/(2B))  -> pulso cuadrado en frecuencia (|f|<=B)
% (Aquí lo usamos ANALÍTICO para que el aliasing se vea perfecto y no "se rompa".)
%
% 1) Tiempo: x(t) continua vs muestreada x[n] vs diezmada y[m]=x[nM]
% 2) Frecuencia analógica (Hz): réplicas de X(f) con fs y con fs/M (aliasing visible si B > fs/(2M))
% 3) Frecuencia digital: superposición |X(e^jw)| vs |Y(e^jw)| en eje nu (ciclos/muestra)

clear; close all; clc;

%% =========================
% PARÁMETROS (TOCA AQUÍ)
%% =========================
fs   = 200;          % Hz (muestreo original)
Ts   = 1/fs;

M    = 2;            % factor de diezmado (entero >=2)
fs_d = fs/M;         % nueva frecuencia de muestreo tras diezmar
Ts_d = 1/fs_d;

B    = 60;           % Hz (semiancho del rect en frecuencia: |f|<=B)
Tobs = 1.5;          % s (ventana temporal para ver sinc)

oversamp = 50;       % oversampling para "continuo": dt = Ts/oversamp
dt = Ts/oversamp;

Ndft = 32768;        % FFT para espectro digital (secuencias)
Krep = 3;            % réplicas a cada lado en Hz (plot analógico)

%% =========================
% 0) Señal "continua" y secuencias
%% =========================
t = -Tobs/2:dt:Tobs/2-dt;

% Señal continua
x = sinc(2*B*t);

% Muestreo a fs
t_s = (-Tobs/2):Ts:(Tobs/2);
x_n = sinc(2*B*t_s);

% Diezmado por M (submuestreo): y[m]=x[nM]
t_d = t_s(1:M:end);
y_m = x_n(1:M:end);

%% =========================
% 1) TIEMPO: continua vs muestreada vs diezmada
%% =========================
figure('Name','(1) Tiempo: continua vs muestreada vs diezmada');
grid on; hold on;

plot(t, x, 'k', 'LineWidth', 1.6);  % continua negra

plot(t_s, x_n, '-o', ...
    'LineWidth', 1.2, ...
    'MarkerSize', 5);

plot(t_d, y_m, '-*', ...
    'LineWidth', 1.6, ...
    'MarkerSize', 8);

xlabel('Tiempo (s)');
ylabel('Amplitud');
title(sprintf('(1) x(t)=sinc(2Bt) | fs=%.1f Hz, M=%d (fs_d=%.1f Hz), B=%.1f Hz', ...
    fs, M, fs_d, B));

legend({'x(t) continua', 'x[n] (fs) -o', 'y[m] (fs/M) -*'}, 'Location','best');

% Zoom para ver bien muestreo vs diezmado
xlim([-(10*Ts_d) (10*Ts_d)]);

%% =========================
% 2) FRECUENCIA ANALÓGICA (Hz): RÉPLICAS ANALÍTICAS (SIEMPRE BIEN, CON ALIASING)
% X(f) = rect(|f|<=B) (normalizada)
% Muestreo ideal: Xs(f) = (1/Ts) * sum_k X(f - k*fs)   (réplicas cada fs)
% Diezmado equivale a muestrear a fs_d: Xsd(f) replicas cada fs_d
%
% Aquí dibujamos SOLO magnitud normalizada (forma), sin preocuparnos por 1/Ts
%% =========================

% Eje de frecuencia analógica para mostrar réplicas
f_plot = linspace(-(Krep+0.5)*fs, (Krep+0.5)*fs, Ndft);

% Rect baseband: 1 si |f|<=B
X0 = double(abs(f_plot) <= B); % pulso cuadrado ideal

% Réplicas con fs (muestreada) y con fs_d (diezmada)
Xs  = zeros(size(f_plot));
Xsd = zeros(size(f_plot));

for k = -Krep:Krep
    Xs  = Xs  + double(abs(f_plot - k*fs)  <= B);
    Xsd = Xsd + double(abs(f_plot - k*fs_d) <= B);
end

% Normaliza para comparar
% Xs  = Xs  / (max(Xs)  + eps);
% Xsd = Xsd / (max(Xsd) + eps);



%% =========================
% (2A) FRECUENCIA ANALÓGICA (Hz): MUETREADA
%% =========================
figure('Name','(2A) Frecuencia analógica (Hz): muestreada');
grid on; hold on;

plot(f_plot, Xs, ...
    'LineWidth', 2.0, ...
    'Color', [0 0.4470 0.7410]);   % azul

xlabel('Frecuencia analógica (Hz)');
ylabel('Magnitud');
title(sprintf('(2A) Réplicas analógicas de la señal muestreada | fs=%.1f Hz, B=%.1f Hz', fs, B));

xlim([-(Krep+0.5)*fs (Krep+0.5)*fs]);

% Guías en múltiplos de fs
yl = ylim;
for k = -Krep:Krep
    xline(k*fs, ':', 'HandleVisibility','off');
end
ylim(yl);


%% =========================
% (2B) FRECUENCIA ANALÓGICA (Hz): DIEZMADA
%% =========================
figure('Name','(2B) Frecuencia analógica (Hz): diezmada');
grid on; hold on;

plot(f_plot, Xsd, ...
    'LineWidth', 1.8, ...
    'Color', [0.8500 0.3250 0.0980]);   % naranja

xlabel('Frecuencia analógica (Hz)');
ylabel('Magnitud');
title(sprintf('(2B) Réplicas analógicas tras diezmado | fs_d=%.1f Hz (M=%d), B=%.1f Hz', fs_d, M, B));

xlim([-(Krep+0.5)*fs (Krep+0.5)*fs]);

% Guías en múltiplos de fs_d
yl = ylim;
for k = -Krep:Krep
    xline(k*fs_d, '--', 'HandleVisibility','off');
end
ylim(yl);


%% =========================
% 3) FRECUENCIA DIGITAL: |X(e^jw)| vs |Y(e^jw)| EN EL MISMO EJE nu
% Aquí sí usamos FFT de las secuencias x[n] y y[m]
%% =========================

Xn = fftshift(fft(x_n, Ndft));
Yn = fftshift(fft(y_m, Ndft));

w  = linspace(-pi, pi, Ndft);
nu = w/(2*pi);   % ciclos/muestra en [-0.5,0.5]

% Xn_mag = abs(Xn); Xn_mag = Xn_mag/(max(Xn_mag)+eps);
% Yn_mag = abs(Yn); Yn_mag = Yn_mag/(max(Yn_mag)+eps);

Xn_mag = abs(Xn);
Yn_mag = abs(Yn); 
figure('Name','(3) Frecuencia digital: muestreada vs diezmada (superpuesto)');
grid on; hold on;

hX = plot(nu, Xn_mag, 'LineWidth', 2.0, 'Color', [0 0.4470 0.7410]);       % azul
hY = plot(nu, Yn_mag, 'LineWidth', 1.8, 'Color', [0.8500 0.3250 0.0980]);  % naranja

xlabel('\nu (ciclos/muestra)');
ylabel('Magnitud (normalizada)');
title(sprintf('(3) Espectro digital superpuesto | fs=%.1f Hz, fs_d=%.1f Hz (M=%d)', fs, fs_d, M));
xlim([-0.5 0.5]);

% Guías: banda que debería quedar para poder diezmar sin alias (en eje de x[n])
xline(-0.5/M, '--', 'HandleVisibility','off');
xline( 0.5/M, '--', 'HandleVisibility','off');

lgd = legend([hX hY], {'Muestreada', 'Diezmada (ensanchada)'}, 'Location','best');
lgd.AutoUpdate = 'off';

fprintf('\n--- Ensanchamiento en eje digital ---\n');
fprintf('En x[n] ocupa aprox B/fs   = %.4f ciclos/muestra\n', B/fs);
fprintf('En y[m] ocupa aprox B/fs_d = %.4f ciclos/muestra = M*(B/fs)\n\n', B/fs_d);

%% =========================
% FIN
%% =========================
