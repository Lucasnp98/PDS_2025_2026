
function ck = DSF_rectificada(V,N)
%recti2(V,N) calcula los 2N+1 coeficientes de Fourier de una señal
%rectificada en doble onda y amplitud V.

if nargin == 0
    V = 1;
    N = 20;
end

k=-N:N; %genera los 2N+1 términos
%Rectificada onda completa
ck=(V/2)*(sinc((k-1)/2)+sinc((k+1)/2)).*exp(-j*pi*k).*cos(k*pi/2);

% Rectificada en media onda
ck=(V/4)*(sinc((k-1)/2)+sinc((k+1)/2)).*exp(-j.*pi.*k./2);