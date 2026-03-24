clear; clc; close all;

%% User-defined constants
N  = 201;      
Ts = 0.1;
T = 0:Ts:N*Ts-0.1;
sigma2 = 4.2;

%% Load data
load('manipulator1.mat');

%% Vector initializzation
requiredVars = {'u1','u1d','y1','u2','u2d','y2'};

%% Acceleration 
u1dd = derivative(u1d,Ts);
u2dd = derivative(u2d,Ts);

%% Plot the inputs and the output of the first and second data set

% --- Figura 1: Primo Set di Dati ---
figure(1);
sgtitle('Dataset 1: Dati Cinematici e Coppia del Manipolatore'); % Titolo generale

subplot(4,1,1);
plot(T', u1, 'b', 'LineWidth', 1.5); 
title('Posizione Angolare (q_m)');
ylabel('Pos [rad]');
grid on;

subplot(4,1,2);
plot(T', u1d, 'b', 'LineWidth', 1.5);
title('Velocità Angolare (dq_m/dt)');
ylabel('Vel [rad/s]');
grid on;

subplot(4,1,3);
plot(T', u1dd, 'b', 'LineWidth', 1.5);
title('Accelerazione Angolare (d^2q_m/dt^2)');
ylabel('Acc [rad/s^2]');
grid on;

subplot(4,1,4);
plot(T', y1, 'r', 'LineWidth', 1.5); % Colore rosso per evidenziare l'output
title('Coppia Applicata (\tau)');
ylabel('Coppia [Nm]');
xlabel('Tempo [s]');
grid on;


% --- Figura 2: Secondo Set di Dati ---
figure(2);
sgtitle('Dataset 2: Dati Cinematici e Coppia del Manipolatore'); % Titolo generale

subplot(4,1,1);
plot(T', u2, 'b', 'LineWidth', 1.5); 
title('Posizione Angolare (q_m)');
ylabel('Pos [rad]');
grid on;

subplot(4,1,2);
plot(T', u2d, 'b', 'LineWidth', 1.5);
title('Velocità Angolare (dq_m/dt)');
ylabel('Vel [rad/s]');
grid on;

subplot(4,1,3);
plot(T', u2dd, 'b', 'LineWidth', 1.5);
title('Accelerazione Angolare (d^2q_m/dt^2)');
ylabel('Acc [rad/s^2]');
grid on;

subplot(4,1,4);
plot(T', y2, 'r', 'LineWidth', 1.5); % Colore rosso per evidenziare l'output
title('Coppia Applicata (\tau)');
ylabel('Coppia [Nm]');
xlabel('Tempo [s]');
grid on;

%% Inverse dynamics 
% Kernel identification
x1 = [u1, u1d, u1dd];
x2 = [u2, u2d, u2dd];

K1i = 10 * Cauchy_kernel(x1, x1, 100); %case 1
K1e = 10 * Cauchy_kernel(x2, x1, 100); %case 1

K2i = 0.1 * Cauchy_kernel(x1, x1, 100); %case 2
K2e = 0.1 * Cauchy_kernel(x2, x1, 100); %case 2

K3i = 10 * Cauchy_kernel(x1, x1, 1); %case 3
K3e = 10 * Cauchy_kernel(x2, x1, 1); %case 3

% Add gaussian noise
g1 = K1e * (K1i + sigma2*eye(N))^-1 * y1;
g2 = K2e * (K2i + sigma2*eye(N))^-1 * y1;
g3 = K3e * (K3i + sigma2*eye(N))^-1 * y1;

%% Plot gmap (Validation)
figure(3);

subplot(3,1,1);
plot(T', y2, 'LineWidth', 1.5); hold on; 
plot(T', g1, 'r--', 'LineWidth', 1.5);
title('Caso 1: \lambda = 10, \beta = 100');
legend('Reale (y_2)', 'Stimato (g_1)');
grid on;

subplot(3,1,2);
plot(T', y2, 'LineWidth', 1.5); hold on; 
plot(T', g2, 'r--', 'LineWidth', 1.5);
title('Caso 2: \lambda = 0.1, \beta = 100');
legend('Reale (y_2)', 'Stimato (g_2)');
grid on;

subplot(3,1,3);
plot(T', y2, 'LineWidth', 1.5); hold on; 
plot(T', g3, 'r--', 'LineWidth', 1.5);
title('Caso 3: \lambda = 10, \beta = 1');
legend('Reale (y_2)', 'Stimato (g_3)');
grid on;

%% Analisi quantitativa dell'accuratezza (Fit Percentage / R^2)
% Media dei dati reali di validazione
y2_mean = mean(y2);

% Calcolo della percentuale di Fit per i 3 casi
% Formula: 100 * (1 - (Norma dell'errore / Norma della variazione totale))
Fit1 = 100 * (norm(y2 - g1) / norm(y2 - y2_mean));
Fit2 = 100 * (norm(y2 - g2) / norm(y2 - y2_mean));
Fit3 = 100 * (norm(y2 - g3) / norm(y2 - y2_mean));

% Stampa i risultati nella Command Window
fprintf('\n--- Accuratezza del Modello (Error Percentage) ---\n');
fprintf('Caso 1 (lambda=10, beta=100): %.2f%%\n', Fit1);
fprintf('Caso 2 (lambda=0.1, beta=100): %.2f%%\n', Fit2);
fprintf('Caso 3 (lambda=10, beta=1): %.2f%%\n', Fit3);

% Calcola l'MSE del tuo modello migliore (Caso 1)
MSE_Caso1 = mean((y2 - g1).^2);
MSE_Caso2 = mean((y2 - g2).^2);
MSE_Caso3 = mean((y2 - g3).^2);

% Confrontalo con la varianza teorica del rumore
fprintf('MSE calcolato: %.4f\n', MSE_Caso1);
fprintf('MSE calcolato: %.4f\n', MSE_Caso2);
fprintf('MSE calcolato: %.4f\n', MSE_Caso3);
fprintf('Varianza teorica del rumore (sigma^2): %.4f\n', sigma2);

%% Feedforward control
% 1. Create reference signal w(t)
w = zeros(N, 1);
w(11:N) = 10; % Gradino da 0 a 10 al tempo 1s

% 2. Calculate derivatives of the reference signal
wd = derivative(w, Ts);
wdd = derivative(wd, Ts);

% 3. Construct input location for feedforward
xf = [w, wd, wdd];

% 4. Compute Feedforward action yf for the 3 cases
% Case 1: lambda = 10, beta = 100
K1f = 10 * Cauchy_kernel(xf, x1, 100);
yf1 = K1f * ((K1i + sigma2*eye(N))^-1) * y1;

% Case 2: lambda = 0.1, beta = 100
K2f = 0.1 * Cauchy_kernel(xf, x1, 100);
yf2 = K2f * ((K2i + sigma2*eye(N))^-1) * y1;

% Case 3: lambda = 10, beta = 1
K3f = 10 * Cauchy_kernel(xf, x1, 1);
yf3 = K3f * ((K3i + sigma2*eye(N))^-1) * y1;

