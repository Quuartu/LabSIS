%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %      
%  Tutorial of the basic functions of the System Identification Toolbox   %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

%% Model generation / Data generation

N = 400; % length of the data
Ts = 1;  % sampling time

% Coefficientes of the actual model (e.g. ARX model):
% numerators and denominators of the model with feedback
Fden = [1 -1.5 0.7];  
Fnum = [0 1 0.5];    
Gden = [1 0 0];
Gnum = [1 0 0];

% Input generation: e.g. sum of sinusoids whose frequencies are in the
% interval [pi*alpha1 pi*alpha2]
alpha1 = 0.1;
alpha2 = 0.2;
k_sin = 2; % number of sinusoids
u = idinput(N,'sine',[alpha1 alpha2],[],k_sin);

% Generate NORMALIZED WGN noise en0 (i.e. en0 has variance equal to one)
en0 = idinput(N,'rgs'); 

% noise variance sigma0^2 of e0=sigma_0*en0
noiseVar = 5;

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


% Estimation of an ARMAX model

% orders of the ARMAX model
% orders_armax(1) = nA
% orders_armax(2) = nB
% orders_armax(3) = nC
% orders_armax(4) = nk
orders_armax = [2 2 2 1];

% ARMAX model estimation
m_armax = armax(data,orders_armax);

% Plot the coefficient of the estimated model
m_armax.a % A(z)
m_armax.b % B(z)
m_armax.c % C(z)
m_armax.NoiseVariance % sigma^2


% Estimation of an OE model

% orders of the OE model
% orders_oe(1) = nB
% orders_oe(2) = nF
% orders_oe(3) = nk
orders_oe = [2 2 1]; 

% OE model estimation
m_oe = oe(data,orders_oe);

% Plot the coefficient of the estimated model
m_oe.b % B(z)
m_oe.f % F(z)
m_oe.NoiseVariance % sigma^2


% Estimation of a Box-Jenkins model

% coefficients of the BJ model
% orders_bj(1) = nB
% orders_bj(2) = nC
% orders_bj(3) = nD
% orders_bj(4) = nF
% orders_bj(5) = nk
orders_bj = [2 2 2 2 1];

% BJ model generation
m_bj = bj(data,orders_bj);

% Plot the coefficient of the estimated model
m_bj.b % B(z)
m_bj.c % C(z)
m_bj.d % D(z) 
m_bj.f % F(z)
m_bj.NoiseVariance % sigma^2


% Estimation of a narx model

% first part of the regression terms: y(t-5)^3 y(t-6)^ 3u(t)^3 u(t-1)^3 
reg1 = polynomialRegressor({'y1','u1'},{5:6,0:1},3);
% first part of the regression temrm: y(t-5)^9 y(t-6)^9 u(t)^9 u(t-1)^9 
reg2 = polynomialRegressor({'y1','u1'},{5:6,0:1},9);
% estimate the narx model
m_narx=nlarx(data,[reg1 reg2]);

%% Bode plot of the estimated input-output transfer functions

bodeplot(m_arx,m_armax,m_oe,m_bj)



%% Kernel-based method

% Estimation of a (approximated) nonparametric BJ model

% orders_arxReg(1) = nA (practical length)
% orders_arxReg(2) = nB (practical length+1)d
% orders_arxReg(3) = nk
orders_arxReg = [50 50 1];

% Define the matrix sigma^2*K^-1
sigma2Kinv=eye(50+50);
Lambda = norm(sigma2Kinv);
R = sigma2Kinv/Lambda;% WARNING: 
% Lambda=||\hat{sigma}^2*K^-1|| and R=\hat{sigma}^2*K^-1/Lambda
% (i.e. Lambda is not the scale factor)
opt = arxOptions;
opt.Regularization.Lambda = Lambda;
opt.Regularization.R = R; 

% Kernel-based estimation 
m_arxReg = arx(data, orders_arxReg, opt);
m_arxReg.a % A(z)
m_arxReg.b % B(z)
m_arxReg.NoiseVariance % sigma^2


 
 