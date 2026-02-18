

%% h[n] causal (solo n>=0)
h = [2 1 0.5];          % h[0]=2, h[1]=1, h[2]=0.5  (causal)

%% entrada x[n] (cualquiera, pero definida para n>=0)
x = [1 3 -2 4 0 2];      % x[0..5]

%% salida por convolución (MATLAB)
y = conv(x,h);

%% Demostración en un n concreto: y[n0] usa SOLO x[n0], x[n0-1], ...
n0 = 4;                  % elige el instante que quieras

% Fórmula (como h es causal): y[n0] = h0*x[n0] + h1*x[n0-1] + h2*x[n0-2]
y_n0 = h(1)*x(n0+1) + h(2)*x(n0) + h(3)*x(n0-1);

fprintf('y[%d] por conv = %.2f\n', n0, y(n0+1));
fprintf('y[%d] por suma = %.2f\n', n0, y_n0);
fprintf('Para y[%d] solo se usan x[%d], x[%d], x[%d]\n', n0, n0, n0-1, n0-2);

%% Dibujos
n_x = 0:length(x)-1;
n_h = 0:length(h)-1;
n_y = 0:length(y)-1;

figure;
subplot(3,1,1); stem(n_x,x,'filled'); grid on; title('x[n]');
subplot(3,1,2); stem(n_h,h,'filled'); grid on; title('h[n] (causal)');
subplot(3,1,3); stem(n_y,y,'filled'); grid on; title('y[n]=x*h');







%% h[n] NO causal (tiene valores para n<0)
% Índices reales: h[-1]=1 , h[0]=2 , h[1]=0.5
h = [1 2 0.5];          % taps
n_h = -1:1;             % índices de h (NO causal)

%% entrada x[n] (cualquiera, definida para n>=0)
x = [1 3 -2 4 0 2];      % x[0..5]
n_x = 0:length(x)-1;

%% salida por convolución
y = conv(x,h);
n_y = n_x(1)+n_h(1) : n_x(end)+n_h(end);

%% Demostración en un n concreto
n0 = 4;

% Fórmula GENERAL: y[n] = sum_k h[k] x[n-k]
% Aquí k = -1,0,1
y_n0 = ...
    h(1)*x(n0-(-1)+1) + ...   % h[-1] * x[n0+1]  --> FUTURO
    h(2)*x(n0-(0)+1)  + ...   % h[0]  * x[n0]
    h(3)*x(n0-(1)+1);         % h[1]  * x[n0-1]

fprintf('y[%d] por conv = %.2f\n', n0, y(n0 - n_y(1) + 1));
fprintf('y[%d] por suma = %.2f\n', n0, y_n0);
fprintf('Para y[%d] se usan x[%d] (futuro), x[%d], x[%d]\n', ...
        n0, n0+1, n0, n0-1);

%% Dibujos
figure;
subplot(3,1,1); stem(n_x,x,'filled'); grid on; title('x[n]');
subplot(3,1,2); stem(n_h,h,'filled'); grid on; title('h[n] (NO causal)');
subplot(3,1,3); stem(n_y,y,'filled'); grid on; title('y[n]=x*h');
