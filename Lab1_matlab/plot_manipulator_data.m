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

%% Plot the inputs and the output of the first and of the second data set
figure(1);
subplot(4,1,1);
plot(T',u1); 
subplot(4,1,2);
plot(T',u1d);
subplot(4,1,3);
plot(T',u1dd);
subplot(4,1,4);
plot(T',y1);

figure(2);
subplot(4,1,1);
plot(T',u2); 
subplot(4,1,2);
plot(T',u2d);
subplot(4,1,3);
plot(T',u2dd);
subplot(4,1,4);
plot(T',y2);

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

%% Plot gmap
figure(3);
subplot(3,1,1);
plot(T',y2); 
hold on 
plot(T',g1);

subplot(3,1,2);
plot(T',y2); 
hold on 
plot(T',g2);

subplot(3,1,3);
plot(T',y2); 
hold on 
plot(T',g3);

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
yf2 = K2f * ((K2i + sigma2*eye(N)^-1) * y1);

% Case 3: lambda = 10, beta = 1
K3f = 10 * Cauchy_kernel(xf, x1, 1);
yf3 = K3f * ((K3i + sigma2*eye(N))^-1) * y1;

% Save yf1 in yf for first Simulink run
yf = yf1; 