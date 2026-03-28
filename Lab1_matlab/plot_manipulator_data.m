clear
clc
close all

%% 1. Initializzation defined constants
N  = 201;      
Ts = 0.1;
T = 0:Ts:N*Ts-0.1;
lambda1=10;
lambda2=0.1;
lambda3=10;
beta1=100;
beta2=100;
beta3=1;
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

% Case 1 lambda=10, Beta=100
KC1i = Cauchy_kernel(x1, x1, beta1);
KC1e = Cauchy_kernel(x2, x1, beta1);
K1i = lambda1 * KC1i; 
K1e = lambda1 * KC1e;
ValueKC1e=mean(diag(KC1e));

% Case 2 lambda=0.1, Beta=100
KC2i = Cauchy_kernel(x1, x1, beta2);
KC2e = Cauchy_kernel(x2, x1, beta2);
K2i = lambda2 * KC2i; 
K2e = lambda2 * KC2e;
ValueKC2e=mean(diag(KC2e));

% Case 3 lambda=10, Beta=1
KC3i = Cauchy_kernel(x1, x1, beta3);
KC3e = Cauchy_kernel(x2, x1, beta3);
K3i = lambda3 * KC3i; 
K3e = lambda3 * KC3e;
ValueKC3e=mean(diag(KC3e));

% Add gaussian noise
g1 = K1e * (K1i + variance*eye(N))^-1 * y1;
g2 = K2e * (K2i + variance*eye(N))^-1 * y1;
g3 = K3e * (K3i + variance*eye(N))^-1 * y1;

%% 6.1 Plot gmap
figure(3)
plot(T', y2, 'b', 'LineWidth', 1.5); hold on; 
plot(T', g1, 'r', 'LineWidth', 1.5);
title('Case 1: \lambda = 10, \beta = 100');
legend('y_2', 'g_1');
grid on;

%% 6.2 Validation
figure(4);

subplot(3,1,1);
plot(T', y2, 'b', 'LineWidth', 1.5); hold on; 
plot(T', g1, 'r', 'LineWidth', 1.5);
title('Case 1: \lambda = 10, \beta = 100');
legend('y_2', 'g_1');
grid on;

subplot(3,1,2);
plot(T', y2,'b', 'LineWidth', 1.5); hold on; 
plot(T', g2, 'r', 'LineWidth', 1.5);
title('Case 2: \lambda = 0.1, \beta = 100');
legend('y_2', 'g_2');
grid on;

subplot(3,1,3);
plot(T', y2,'b', 'LineWidth', 1.5); hold on; 
plot(T', g3, 'r', 'LineWidth', 1.5);
title('Case 3: \lambda = 10, \beta = 1');
legend('y_2', 'g_3');
grid on;

%% 7. Feedforward control
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
K1f = lambda1 * Cauchy_kernel(xf, x1, beta1);
yf1 = K1f * ((K1i + variance*eye(N))^-1) * y1;

% Case 2
K2f = lambda2 * Cauchy_kernel(xf, x1, beta2);
yf2 = K2f * ((K2i + variance*eye(N))^-1) * y1;

% Case 3
K3f = lambda3 * Cauchy_kernel(xf, x1, beta3);
yf3 = K3f * ((K3i + variance*eye(N))^-1) * y1;

%% 8. Automatic runs in Simulink (switch: 0, yf1, yf2, yf3)
% Minimal flow only: set case, simulate, plot figures 5-8.

modelName = 'scheme';
switchBlock = [modelName '/Switch'];
feedforwardBlock = [modelName '/Feedforward'];
angularScopeBlock = [modelName '/Angular position'];

% From Workspace expects [t yfi], ensure column vectors with same time base.
t = T(:);
yf1 = yf1(:);
yf2 = yf2(:);
yf3 = yf3(:);

load_system(modelName);
set_param(angularScopeBlock, 'DataLogging', 'on');

scopeLabels = {'switch_0', 'switch_yf1', 'switch_yf2', 'switch_yf3'};
switchStates = {'1', '0', '0', '0'};
ffVars = {'yf1', 'yf1', 'yf2', 'yf3'};
plotTitles = {'Without FeedbackForward', 'yf1', 'yf2', 'yf3'};

for k = 1:4
	set_param(feedforwardBlock, 'VariableName', sprintf('[t %s]', ffVars{k}));
	set_param(switchBlock, 'sw', switchStates{k});
	scopeVarName = ['ScopeData_' scopeLabels{k}];
	set_param(angularScopeBlock, 'DataLoggingVariableName', scopeVarName);

	simOut = sim(modelName, 'ReturnWorkspaceOutputs', 'on');

	scopeData = simOut.get(scopeVarName);
	if isa(scopeData, 'Simulink.SimulationData.Dataset')
		tPlot = scopeData{1}.Values.Time(:);
		yPlot = squeeze(scopeData{1}.Values.Data);
		if isvector(yPlot)
			yPlot = yPlot(:);
		end
		if scopeData.numElements >= 2
			y2 = squeeze(scopeData{2}.Values.Data);
			if isvector(y2)
				y2 = y2(:);
			end
			yPlot = [yPlot y2];
		end
	else
		tPlot = scopeData.time(:);
		yPlot = scopeData.signals.values;
	end

	figure(k + 4);
	clf;
	plot(tPlot, yPlot(:,1), 'b', 'LineWidth', 1.5); hold on;
	plot(tPlot, yPlot(:,2), 'r--', 'LineWidth', 1.5);
	legend('Single-link manipulator', 'Reference signal', 'Location', 'best');
	grid on;
	xlabel('Time [s]');
	ylabel('Position [rad]');
	title(plotTitles{k});
end
