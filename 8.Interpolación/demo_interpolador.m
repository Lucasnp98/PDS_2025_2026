%% Interpolación por inserción de ceros + filtro (demo visual con toggle de ganancia L)
clear; close all; clc;

%% Parámetros
B  = 30;        % "B" en sinc(B*t)
Fs = 100;       % Hz (muestreo original)
L  = 2;         % factor de interpolación (inserta L-1 ceros)

usar_ganancia_L = true;   % <---- CAMBIA a false para ver el efecto

Ts  = 1/Fs;
Fs_i = Fs*L;
Ts_i = 1/Fs_i;

%% 1) Secuencia muestreada x[n] = sinc(B*t)
Tobs = 0.6;
t_s  = -Tobs/2 : Ts : Tobs/2;
x_n  = sinc(B*t_s);

%% 2) Inserción de ceros: u[k]
Nu  = L * length(x_n);
u_k = zeros(1, Nu);
u_k(1:L:end) = x_n;

t_i = t_s(1) + (0:Nu-1)*Ts_i;

%% 3) Interpolación (LPF ideal aprox) para "rellenar" los ceros
nu_c = 1/(2*L);

Lfilt = 401;
if mod(Lfilt,2)==0, Lfilt = Lfilt+1; end
n = -(Lfilt-1)/2:(Lfilt-1)/2;

h_ideal = 2*nu_c * sinc(2*nu_c*n);
w = hann(Lfilt)';
h = h_ideal .* w;

% Normaliza a DC=1
h = h / sum(h);

% Toggle de ganancia L (DC pasa de 1 -> L)
if usar_ganancia_L
    h = L * h;
end

y_k = conv(u_k, h, 'same');

%% ===== (A) Plot en TIEMPO (3 curvas superpuestas) =====
figure; hold on; grid on;

plot(t_s, x_n, '-o', 'LineWidth', 1.2, 'MarkerSize', 5);
plot(t_i, u_k, '-*', 'LineWidth', 1.0, 'MarkerSize', 5);
plot(t_i, y_k, '-x', 'LineWidth', 1.6, 'MarkerSize', 5);

xlabel('Tiempo (s)');
ylabel('Amplitud');
ttl = sprintf('Interpolación L=%d | usar\\_ganancia\\_L=%d', L, usar_ganancia_L);
title(ttl);
legend({'x[n] (Fs)  -o', 'u[k] (ceros)  -*', 'y[k] (interpolada)  -x'}, 'Location','best');
xlim([-(10*Ts) (10*Ts)]);

% Check numérico: ¿cuánto coincide y[k] en los instantes originales k=nL?
idx = 1:L:length(y_k);
gain_est = mean(y_k(idx) ./ (x_n + eps));
fprintf('Ganancia media en instantes k=nL: %.4f (esperado ~1 si usar_ganancia_L=true, ~1/L si false)\n', gain_est);

%% ===== (B) Espectro digital (para ver imágenes y escala correcta) =====
Ndft = 32768;

Xn = fftshift(fft(x_n, Ndft));
U  = fftshift(fft(u_k, Ndft));
Y  = fftshift(fft(y_k, Ndft));

nu = (-Ndft/2 : Ndft/2-1) / Ndft;

% Escalado para comparar amplitudes (importante!)
Xn_mag = abs(Xn) / length(x_n);
U_mag  = abs(U)  / length(u_k);
Y_mag  = abs(Y)  / length(y_k);

figure; hold on; grid on;
plot(nu, Xn_mag, 'LineWidth', 1.6);
plot(nu, U_mag,  'LineWidth', 1.2);
plot(nu, Y_mag,  'LineWidth', 1.6);

xlabel('\nu (ciclos/muestra)');
ylabel('Magnitud (FFT escalada)');
title(sprintf('Espectro digital (escalado) | L=%d | usar\\_ganancia\\_L=%d', L, usar_ganancia_L));
legend('|X_n(\nu)|','|U(\nu)| (imágenes)','|Y(\nu)| (LPF)','Location','best');
xlim([-0.5 0.5]);

% Guías en k/L
for k = -ceil(L/2):ceil(L/2)
    xline(k/L, ':', 'HandleVisibility','off');
end

% Corte ideal del interpolador
xline( 1/(2*L), '--', 'HandleVisibility','off');
xline(-1/(2*L), '--', 'HandleVisibility','off');
