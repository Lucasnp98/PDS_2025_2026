%% =========================================================
function X = dtft_vis(x, n, w)
    % x: vector de muestras de la secuencia
    % n: índices correspondientes a x (mismo tamaño que x)
    % w: vector de frecuencias (rad/muestra)
    X = zeros(size(w));
    for k = 1:numel(w)
        X(k) = sum( x .* exp(-1j*w(k)*n) );
    end
end
