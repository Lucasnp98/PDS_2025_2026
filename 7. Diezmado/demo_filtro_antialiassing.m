%% Sinc^2 y su FT en frecuencia digital + antialiasing + diezmado
clear; close all; clc;

B  = 30;      % Hz
Fs = 100;     % Hz
Ts = 1/Fs;

M = 2;                    % factor de diezmado
filtro_antialiassing = true;

%% --- Longitud de señal ---
T  = 12;                      % segundos
N  = round(T*Fs);
if mod(N,2)==1, N = N+1; end

t = (-N/2 : N/2-1) * Ts;

%% --- Señal sinc^2 (BW = [-B, B] Hz) ---
x = (B*sinc(B*t)).^2;

%% --- Antialiasing LPF ideal (frecuencia digital) ---
fc_d = 1/(2*M);               % corte en frecuencia digital (f/Fs), ciclos/muestra

L = 801;                      % longitud FIR (impar)
if mod(L,2)==0, L=L+1; end
n = -(L-1)/2:(L-1)/2;

h_ideal = 2*fc_d * sinc(2*fc_d*n);   % LPF ideal (en f_d)
w = hann(L)';                         % ventana
h = h_ideal .* w;

h = h / sum(h);                       % normalización ganancia DC=1

%% --- Aplicar antialiasing ANTES de cualquier FFT / diezmado (si procede) ---
x_f = x;                               % por claridad
if filtro_antialiassing
    x_f = conv(x, h, 'same');          % (didáctico, con bordes)
    % Alternativa más "limpia" (no causal): x_f = filtfilt(h, 1, x);
end

%% --- FFT de la señal (original y/o filtrada) ---
X  = fftshift(fft(x));                 % espectro de la señal original
Xmag = abs(X);  Xmag = Xmag / max(Xmag);

Xf = fftshift(fft(x_f));               % espectro de la señal que se va a diezmar
Xfmag = abs(Xf); Xfmag = Xfmag / max(Xfmag);

f_d_sig = (-N/2 : N/2-1) / N;          % frecuencia digital (ciclos/muestra)

%% --- Respuesta en frecuencia del filtro (frecuencia digital) ---
Nfft = 4096;
H = fftshift(fft(h, Nfft));
Hmag = abs(H); Hmag = Hmag / max(Hmag);

f_d_filt = (-Nfft/2 : Nfft/2-1) / Nfft;

%% --- Plot 1: Señal y filtro en frecuencia digital ---
figure; hold on; grid on;

plot(f_d_sig,  Xmag,  'LineWidth', 1.4);     % señal original
plot(f_d_filt, Hmag,  'LineWidth', 1.4);     % filtro
if filtro_antialiassing
    plot(f_d_sig, Xfmag, 'LineWidth', 1.4);  % señal filtrada (opcional)
end

xlabel('f_d = f / F_s (ciclos/muestra)');
ylabel('Magnitud normalizada');
title('Señal sinc^2 y filtro antialiasing en frecuencia digital');

if filtro_antialiassing
    legend('|X(f_d)| (original)','|H(f_d)|','|X_f(f_d)| (filtrada)','Location','Best');
else
    legend('|X(f_d)| (original)','|H(f_d)|','Location','Best');
end

% Banda de la señal (original) en f_d
f_dc_sig = B/Fs;
xline(+f_dc_sig,'--','Signal BW','LabelVerticalAlignment','bottom');
xline(-f_dc_sig,'--');

% Corte del filtro en f_d
f_dc_filt = 1/(2*M);
xline(+f_dc_filt,'--','Filter cut','LabelVerticalAlignment','bottom');
xline(-f_dc_filt,'--');

xlim([-0.5 0.5]);

%% --- Diezmado + espectro en frecuencia digital referida a Fs_new ---
Fs_new = Fs / M;

x_dec = x_f(1:M:end);                    % diezmado de la señal (filtrada o no)
Ndec  = length(x_dec);

Xdec = fftshift(fft(x_dec));
Xdec_mag = abs(Xdec);
Xdec_mag = Xdec_mag / max(Xdec_mag);

f_d_dec = (-Ndec/2 : Ndec/2-1) / Ndec;   % f_d = f / Fs_new (ciclos/muestra)

figure; grid on;
plot(f_d_dec, Xdec_mag, 'LineWidth', 1.3);
xlabel('f_d = f / F_{s,new} (ciclos/muestra)');
ylabel('Magnitud normalizada');
if filtro_antialiassing
    title('Espectro tras filtrar (LPF) y diezmar (frecuencia digital recalculada)');
else
    title('Espectro tras diezmar (sin antialiasing) (frecuencia digital recalculada)');
end
xlim([-0.5 0.5]);
