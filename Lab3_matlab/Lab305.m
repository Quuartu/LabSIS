clc
clear all
close all

%% DATA PLOTTING
load('data.mat');
figure(1);

plot(u);
hold on;
plot(y);
title("Input and Output plot");
legend('Input (u)', 'Output (y)');
grid on;

%% Kernel Matrix
nB = 40;

% Parameter to change for testing cases
lambda = 2;
beta = 0.75;

K_base = zeros(nB, nB);

% --- TC KERNEL (Tuned Correlated) ---
for i = 1:nB
    for j = 1:nB
        K_base(i,j) = beta^max(i, j);
    end
end

% --- DI KERNEL (Diagonal) ---
% for i = 1:nB
%     K_base(i,i) = beta^i;     
% end

% --- SS KERNEL (Stable Spline) ---
% c = beta^(1/3);
% for i = 1:nB
%    for j = 1:nB 
%        K_base(i,j) = ((c^(i+j) * min(c^i, c^j)) / 2) - (min(beta^i, beta^j) / 6);
%    end
% end

K = lambda * K_base;
L = chol(K, 'lower'); 

theta_realizations = zeros(nB, 10);

for i = 1:10
    v = randn(nB,1);
    theta_realizations(:, i) = L * v;
end

figure(2);
plot(theta_realizations);
title("10 \theta Realizations")

xlabel('k');
ylabel('Coefficient Amplitude');
grid on;

num_realizations = size(theta_realizations, 2);
legend_labels = cell(1, num_realizations);

for i = 1:num_realizations
   legend_labels{i} = sprintf('Realization %d', i);
end

legend(legend_labels, 'Location', 'eastoutside');

%% Kernel Based Estimation

data = iddata(y, u, 1); 

% orders_arxReg(1) = nA (practical length)
% orders_arxReg(2) = nB (practical length+1)d
% orders_arxReg(3) = nk
orders_arxReg = [0 nB 1]; % delay (t-1)

% Define the matrix sigma^2*K^-1
sigma2Kinv = (0.58^2) * K^-1; 
%sigma2Kinv = eye(1+40);
Lambda = norm(sigma2Kinv);
R = sigma2Kinv/Lambda;% WARNING: 
% Lambda=||\hat{sigma}^2*K^-1|| and R=\hat{sigma}^2*K^-1/Lambda
% (i.e. Lambda is not the scale factor)
opt = arxOptions;
opt.Regularization.Lambda = Lambda;
opt.Regularization.R = R; 

% Kernel-based estimation 
m_arxReg = arx(data, orders_arxReg, opt);
m_arxReg.a; % A(z)
m_arxReg.b; % B(z)
m_arxReg.NoiseVariance; % sigma^2

figure(3);
plot(theta0);
hold on;
plot(m_arxReg.b(2:end), 'r--', 'LineWidth', 1.5);

title('Comparison: True vs. Estimated Impulse Response');
xlabel('k');
ylabel('Amplitude');

legend('\theta_0 (True)', '\theta (Estimated)', 'Location', 'best');
grid on;

%% Confidence Interval
x_1 = [0; u(1:length(u)-1)];
x_2 = zeros(1,nB);

Phi = toeplitz(x_1,x_2);

P_hat = K-K*Phi'*(Phi*K*Phi'+m_arxReg.NoiseVariance*eye(100))^-1 *Phi*K;

std_dev = sqrt(diag(P_hat));

b_hat = m_arxReg.b(2:end)'; 

upper_bound = b_hat + 1.96 * std_dev;
lower_bound = b_hat - 1.96 * std_dev;

figure(4);
plot(theta0, 'b', 'LineWidth', 1.5); 
hold on;
plot(b_hat, 'g--', 'LineWidth', 1.5);
plot(upper_bound, 'r--', 'LineWidth', 0.5);
plot(lower_bound, 'r--', 'LineWidth', 0.5);
title('Impulse response & confidence intervals');
legend('\theta_0', '\theta', 'Upper Interval', 'Lower Interval');
grid on;