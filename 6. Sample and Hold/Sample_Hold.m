%% ZOH (Sample & Hold) en frecuencia - 3 graficos "bien"
% Añadido: booleano para representar en dB o en magnitud lineal (unidades naturales)

clear; close all; clc;

%% =======================
% Parametros
%% =======================
use_dB = true;      % <-- CAMBIA a false para magnitud lineal

fs   = 80;          % Hz
Ts   = 1/fs;

B    = 8;           % Hz 
Krep = 2;           % <-- menos replicas suele verse mas claro

Tobs = 6.0;         % s
dt   = 1e-4;        % s

t = -Tobs/2:dt:Tobs/2-dt;

%% =======================
% Señal en tiempo: sinc^2  -> espectro triangular (magnitud)
%% =======================
x = (sinc(B*t)).^2;

%% =======================
% CTFT aproximada via FFT (escala dt)
%% =======================
N    = numel(t);
Nfft = 2^nextpow2(4*N);
f    = (-Nfft/2:Nfft/2-1) * (1/(dt*Nfft));  % Hz

X = fftshift(fft(x, Nfft)) * dt;

%% =======================
% (1) Espectro analogico (sin replicas)
%% =======================
magX = abs(X);

% Normalizacion SOLO para grafico 1: a DC (para ver la forma)
[~, i0] = min(abs(f));
magX1 = magX / (magX(i0) + eps);

%% =======================
% Espectro muestreado ideal con replicas: Xs(f) = (1/Ts)*sum_k X(f-kfs)
%% =======================
Xs = zeros(size(f));
for k = -Krep:Krep
    X_shift = interp1(f, X, f - k*fs, 'linear', 0); % complejo
    Xs = Xs + X_shift;
end
Xs = (1/Ts) * Xs;
magXs = abs(Xs);

%% =======================
% ZOH: H(f)=Ts*sinc(fTs)*exp(-j*pi*fTs) -> |H|=Ts*|sinc(fTs)|
%% =======================
Hzoh = Ts * sinc(f*Ts) .* exp(-1j*pi*f*Ts);
Hzoh_mag = abs(Hzoh);

%% =======================
% Overlay de sinc para grafico 2 (SIN deformar): escala solo vertical
%% =======================
idxBB = (f >= -fs/2) & (f <= fs/2);
peakBB = max(magXs(idxBB));
Hzoh_overlay = Hzoh_mag * (peakBB / max(Hzoh_mag + eps));

%% =======================
% Helpers de representacion (dB vs lineal)
%% =======================
toPlot = @(x) x;
yLabel = 'Magnitud';

if use_dB
    toPlot = @(x) 20*log10(x + 1e-15);
    yLabel = 'Magnitud (dB)';
end

%% =======================
% Rango de frecuencia mostrado
%% =======================
fmax = (Krep+1)*fs;

%% =======================
% (1) Grafico: |X(f)| analogo
%% =======================
figure('Name','1) |X(f)| analogo (sin replicas)'); grid on; hold on;
plot(f, toPlot(magX1), 'LineWidth', 1.8);
xlim([-3*B 3*B]);

xlabel('Frecuencia (Hz)');
ylabel(yLabel);
title('1) Espectro analógico: |X(f)| de x(t)=sinc^2(Bt) (sin réplicas)');

if use_dB, ylim([-120 10]); else, ylim([0 1.05]); end

%% =======================
% (2) Grafico: replicas + sinc overlay
%% =======================
figure('Name','2) Réplicas + sinc (ZOH) superpuesta'); grid on; hold on;
plot(f, toPlot(magXs), 'LineWidth', 1.5);
plot(f, toPlot(Hzoh_overlay), '--', 'LineWidth', 1.8);

xlim([-fmax fmax]);
xlabel('Frecuencia (Hz)');
ylabel(yLabel);
title('2) |X_s(f)| con réplicas en k f_s + envolvente sinc del ZOH (overlay)');

legend('|X_s(f)| (réplicas)', '|H_{ZOH}(f)| overlay', 'Location','best');

for m = -Krep:Krep
    xline(m*fs, ':');
end

if use_dB, ylim([-120 20]); else, ylim([0 max([magXs,Hzoh_overlay],[],'all')*1.05]); end

%% =======================
% (3) Grafico: resultado Xzoh = Xs * Hzoh (distorsion)
%% =======================
Xzoh = Xs .* Hzoh;
magXzoh = abs(Xzoh);

figure('Name','3) Resultado: réplicas * sinc (distorsión ZOH)'); grid on; hold on;
plot(f, toPlot(magXzoh), 'LineWidth', 1.8);

xlim([-fmax fmax]);
xlabel('Frecuencia (Hz)');
ylabel(yLabel);
title('3) |X_{ZOH}(f)| = |X_s(f) H_{ZOH}(f)| (distorsión por Sample & Hold)');

for m = -Krep:Krep
    xline(m*fs, ':');
end

if use_dB
    ylim([-120 20]);
else
    ylim([0 max(magXzoh,[],'all')*1.05]);
end
