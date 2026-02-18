% =========================================================
% DTFT de x[n] = a^n u[n]  |a|<1
% Variando el parametro a
% =========================================================

clear; clc; close all;

%% Eje de frecuencia digital
w = linspace(-pi, pi, 2000);   % omega en [-pi, pi]

%% Valores de a a comparar
a_values = [0.2 0.5 0.8 0.95];

figure('Name','DTFT de a^n u[n] para distintos a','NumberTitle','off');

for k = 1:length(a_values)

    a = a_values(k);

    % DTFT analitica:
    % X(e^{jw}) = 1 / (1 - a e^{-jw})
    X = 1 ./ (1 - a*exp(-1j*w));

    % Modulo
    subplot(length(a_values),1,k)
    plot(w, abs(X), 'LineWidth', 1.5);
    grid on;

    title(sprintf('a = %.2f', a));
    xlabel('\omega (rad/muestra)');
    ylabel('|X(e^{j\omega})|');

end
