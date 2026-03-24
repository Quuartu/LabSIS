function u_dot = derivative(u, Ts)


% filter coefficients

a = [1 -0.5];
b= 1/sum(abs(dimpulse(tf([1 0],[1 -0.5],-1),100)));
% computing the derivative after filtering the signal
u_dot = diff(filter(b,a,u))./Ts;
u_dot = [0; u_dot];
end


