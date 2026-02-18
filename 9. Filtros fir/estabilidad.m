% =========================================================
% SECUENCIAS EXPONENCIALES: GRAFICA + ESTABILIDAD (CON BUCLE)
% =========================================================

clear; clc; close all;

%% Eje temporal
n = -10:10;

%% Escalones


%% Definicion de secuencias (handles)
seqs = {
    @(n) (0.5).^n   .* (n>=0),      '(1/2)^n u[n]'
    @(n) (0.5).^(-n).* (n>=0),      '(1/2)^{-n} u[n]'
    @(n) (2).^n     .* (n>=0),      '2^n u[n]'
    @(n) (2).^(-n)  .* (n>=0),      '2^{-n} u[n]'
    @(n) (2).^(-n)  .* (n<=0),      '2^{-n} u[-n]'
};

%% Bucle: grafica y suma
disp('-------------------------------------------');
disp('SUMA ABSOLUTA (aprox)  ->  ESTABILIDAD');
disp('-------------------------------------------');

for k = 1:size(seqs,1)

    f     = seqs{k,1};
    label = seqs{k,2};

    x = f(n);

    S = sum(abs(x));

    fprintf('%-20s  sum = %12.3e\n', label, S);

    figure;
    stem(n, x, 'filled');
    grid on;
    title(label);
    xlabel('n'); ylabel('x[n]');

end
