
function cn = DSF_pulso(D,V,N)

if nargin == 0
    D = .35;
    V = 1;
    N = 20;
end

k = -N:N;

cn = V*D*exp(-j*pi*k*D).*sinc(k*D);

if nargin == 0
    figure
    stem(k,cn);
end