function DSF
    V = 1;
    N = 10;          % mitad del número de coeficientes (±N)
    D = 0.5;         % duty cycle
    T = 5;           % periodo fundamental
    w0 = 2*pi/T;
    t = linspace(-5*T, 5*T, 1e4);

    % Coeficientes c_k para k = -N:N (suponiendo que DSF_pulso los devuelve en ese orden)
    ck = DSF_pulso(D, V, N);
    ck = ck(:).';               % asegúrate de que es fila [1 x (2N+1)]

    k = -N:N;                   % índices armónicos

    x = 0;                      % <-- como pides: escalar, no vector
    for j = -N:N
        % OJO: producto, no suma; y operaciones elemento a elemento con t (.*)
        x = x + ck(j + N + 1) .* exp(1i*j*w0 .* t);
    end

    % Si esperas señal real (pulso real), toma la parte real al final
    x = real(x);

    figure;
    plot(t, x, 'LineWidth', 1.2); grid on
    xlabel('Tiempo [s]');
    ylabel('Amplitud');
    title('Síntesis de una señal pulso');
end