clear all
close all
clc

%% Model generation / Data generation

N = 200; % length of the data
Ts = 1;  % sampling time

% Coefficientes of the actual model (e.g. ARX model):
% numerators and denominators of the model with feedback
Fden = [1 -0.96  0.97];  
Fnum = [0 2.99 -0.2];    
Gden = [1 -0.96 0.97];
Gnum = [1 0 0];

% Input generation: e.g. sum of sinusoids whose frequencies are in the
% interval [pi*alpha1 pi*alpha2]
alpha1 = 1/8;
alpha2 = 1/2;
k_sin = 5; % number of sinusoids
u = idinput(N,'sine',[alpha1 alpha2],[-4 4],k_sin);

% Generate NORMALIZED WGN noise en0 (i.e. en0 has variance equal to one)
en0 = idinput(N,'rgs'); 

% noise variance sigma0^2 of e0=sigma_0*en0
noiseVar = 4.6^2;

% Creation of the actual model (Ts is the sampling time, this parameter
% is optional and the default value is Ts = 1)
m0 = idpoly([],Fnum,Gnum,Gden,Fden,noiseVar,Ts); 

% Creation of the iddata object (a container for all the identification 
% data), in this case we store only the input data u
u = iddata([],u,Ts); 

% Creation of the iddata object, in this case we store only the normalized 
% noise data en0
en0 = iddata([],en0,Ts); 

% Generation of the output data given model, input and noise
y = sim(m0, [u en0]);

% Creation of an iddata object for the output and input data 
data = iddata(y,u);

%% PEM method

% Estimation of an ARX model

% orders of the ARX model
% orders_arx(1) = nA
% orders_arx(2) = nB
% orders_arx(3) = nk
orders_arx = [2 2 1]; 

% ARX model estimation
m_arx = arx(data,orders_arx);

% Plot the coefficients of the estimated model
m_arx.a % A(z)
m_arx.b % B(z)
m_arx.NoiseVariance % sigma^2

figure(1)
bodeplot(m0, m_arx);
legend('True M0', 'M1 (ARX)', 'Location', 'best');

% Estimation of an ARMAX model

% orders of the ARMAX model
% orders_armax(1) = nA
% orders_armax(2) = nB
% orders_armax(3) = nC
% orders_armax(4) = nk
orders_armax = [2 2 1 1];

% ARMAX model estimation
m_armax = armax(data,orders_armax);

% Plot the coefficient of the estimated model
m_armax.a % A(z)
m_armax.b % B(z)
m_armax.c % C(z)
m_armax.NoiseVariance % sigma^2

figure(2)
bodeplot(m0, m_armax);
legend('True M0', 'M2 (ARMAX)', 'Location', 'best');

% Estimation of an OE model

% orders of the OE model
% orders_oe(1) = nB
% orders_oe(2) = nF
% orders_oe(3) = nk
orders_oe = [2 1 1]; 

% OE model estimation
m_oe = oe(data,orders_oe);

% Plot the coefficient of the estimated model
m_oe.b % B(z)
m_oe.f % F(z)
m_oe.NoiseVariance % sigma^2

figure(3)
bodeplot(m0, m_oe);
legend('True M0', 'M3 (OE)', 'Location', 'best');

% Estimation of a Box-Jenkins model

% coefficients of the BJ model
% orders_bj(1) = nB
% orders_bj(2) = nC
% orders_bj(3) = nD
% orders_bj(4) = nF
% orders_bj(5) = nk
orders_bj = [2 1 1 2 1];

% BJ model generation
m_bj = bj(data,orders_bj);

% Plot the coefficient of the estimated model
m_bj.b % B(z)
m_bj.c % C(z)
m_bj.d % D(z) 
m_bj.f % F(z)
m_bj.NoiseVariance % sigma^2

figure(4)
bodeplot(m0, m_bj);
legend('True M0', 'M4 (BJ)', 'Location', 'best');

%% Confidence intervals

theta0 = [-0.96; 0.97; 2.99; -0.2]; 

theta_est = [m_arx.a(2); m_arx.a(3); m_arx.b(2); m_arx.b(3)];
sigma_est2 = m_arx.NoiseVariance;

y_val = y.OutputData;
u_val = u.InputData;

Psi = zeros(4, N);
for t = 3:N
    Psi(:, t) = [-y_val(t-1); -y_val(t-2); u_val(t-1); u_val(t-2)];
end

R_N = (Psi * Psi') / N;
P_hat_N = sigma_est2 * inv(R_N);

for i = 1:4
    cert_N = 1.96 * sqrt(P_hat_N(i,i) / N);
    err_N = abs(theta_est(i) - theta0(i));
    in_N = err_N <= cert_N;

    fprintf('  [N = 200]  Error: %.4f | Confidence: %.4f | Inclusion: %d\n', err_N, cert_N, in_N);
end

