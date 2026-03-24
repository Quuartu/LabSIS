clear; clc; close all;

%% 1. Initializzation defined constants
N  = 201;      
Ts = 0.1;
T = 0:Ts:N*Ts-0.1;
variance = 4.2;

%% 2. Data Loading
load('manipulator1.mat');

%% Vector initializzation
requiredVars = {'u1','u1d','y1','u2','u2d','y2'};

%% 3. Acceleration Computation
u1dd = derivative(u1d,Ts);
u2dd = derivative(u2d,Ts);

%% 4. Data Plotting 

figure(1);
sgtitle('Dataset 1'); 

subplot(4,1,1);
plot(T', u1, 'b', 'LineWidth', 1.5); 
title('q_m');
ylabel('Pos [rad]');
grid on;

subplot(4,1,2);
plot(T', u1d, 'b', 'LineWidth', 1.5);
title('dq_m/dt');
ylabel('Vel [rad/s]');
grid on;

subplot(4,1,3);
plot(T', u1dd, 'b', 'LineWidth', 1.5);
title('d^2q_m/dt^2');
ylabel('Acc [rad/s^2]');
grid on;

subplot(4,1,4);
plot(T', y1, 'r', 'LineWidth', 1.5);
title('\tau');
ylabel('Torque [Nm]');
xlabel('Time [s]');
grid on;


figure(2);
sgtitle('Dataset 2'); 

subplot(4,1,1);
plot(T', u2, 'b', 'LineWidth', 1.5); 
title('q_m');
ylabel('Pos [rad]');
grid on;

subplot(4,1,2);
plot(T', u2d, 'b', 'LineWidth', 1.5);
title('dq_m/dt');
ylabel('Vel [rad/s]');
grid on;

subplot(4,1,3);
plot(T', u2dd, 'b', 'LineWidth', 1.5);
title('d^2q_m/dt^2');
ylabel('Acc [rad/s^2]');
grid on;

subplot(4,1,4);
plot(T', y2, 'r', 'LineWidth', 1.5); 
title('\tau');
ylabel('Torque [Nm]');
xlabel('Time [s]');
grid on;

%% 5. Inverse dynamics 
% Kernel identification
x1 = [u1, u1d, u1dd];
x2 = [u2, u2d, u2dd];

% Case 1
K1i = 10 * Cauchy_kernel(x1, x1, 100); 
K1e = 10 * Cauchy_kernel(x2, x1, 100); 

% Case 2
K2i = 0.1 * Cauchy_kernel(x1, x1, 100); 
K2e = 0.1 * Cauchy_kernel(x2, x1, 100); 

% Case 3
K3i = 10 * Cauchy_kernel(x1, x1, 1); 
K3e = 10 * Cauchy_kernel(x2, x1, 1); 

% Add gaussian noise
g1 = K1e * (K1i + variance*eye(N))^-1 * y1;
g2 = K2e * (K2i + variance*eye(N))^-1 * y1;
g3 = K3e * (K3i + variance*eye(N))^-1 * y1;

%% 6. Plot gmap (Validation)
figure(3);

subplot(3,1,1);
plot(T', y2, 'LineWidth', 1.5); hold on; 
plot(T', g1, 'r--', 'LineWidth', 1.5);
title('Case 1: \lambda = 10, \beta = 100');
legend('y_2', 'g_1');
grid on;

subplot(3,1,2);
plot(T', y2, 'LineWidth', 1.5); hold on; 
plot(T', g2, 'r--', 'LineWidth', 1.5);
title('Case 2: \lambda = 0.1, \beta = 100');
legend('y_2', 'g_2');
grid on;

subplot(3,1,3);
plot(T', y2, 'LineWidth', 1.5); hold on; 
plot(T', g3, 'r--', 'LineWidth', 1.5);
title('Case 3: \lambda = 10, \beta = 1');
legend('y_2', 'g_3');
grid on;

%% 7. Quantitative Analysis of Error
% Mean of y2
y2_mean = mean(y2);

% Formula: 100 * (1 - (Norm Error / Norm of variance))
Err1 = 100 * (norm(y2 - g1) / norm(y2 - y2_mean));
Err2 = 100 * (norm(y2 - g2) / norm(y2 - y2_mean));
Err3 = 100 * (norm(y2 - g3) / norm(y2 - y2_mean));

fprintf('\n--- Error Percentage ---\n');
fprintf('Case 1: %.2f%%\n', Err1);
fprintf('Case 2: %.2f%%\n', Err2);
fprintf('Case 3: %.2f%%\n', Err3);

% Compute MSE
MSE_Case1 = mean((y2 - g1).^2);
MSE_Case2 = mean((y2 - g2).^2);
MSE_Case3 = mean((y2 - g3).^2);

fprintf('\n--- MSE ---\n');
fprintf('Case 1: %.4f\n', MSE_Case1);
fprintf('Case 2: %.4f\n', MSE_Case2);
fprintf('Case 3: %.4f\n', MSE_Case3);

%% 8. Feedforward control
% Creation of signal w(t): step from 0 to 10
w = zeros(N, 1);
w(11:N) = 10; 

% Calculate derivatives of the reference signal
wd = derivative(w, Ts);
wdd = derivative(wd, Ts);

% Construct input location for feedforward
xf = [w, wd, wdd];

% Compute Feedforward action yf for the 3 cases
% Case 1
K1f = 10 * Cauchy_kernel(xf, x1, 100);
yf1 = K1f * ((K1i + variance*eye(N))^-1) * y1;

% Case 2
K2f = 0.1 * Cauchy_kernel(xf, x1, 100);
yf2 = K2f * ((K2i + variance*eye(N))^-1) * y1;

% Case 3
K3f = 10 * Cauchy_kernel(xf, x1, 1);
yf3 = K3f * ((K3i + variance*eye(N))^-1) * y1;

% Set in scheme.slx: yf1 for the first run. Change it to yf2 or yf3 to test the other cases.

