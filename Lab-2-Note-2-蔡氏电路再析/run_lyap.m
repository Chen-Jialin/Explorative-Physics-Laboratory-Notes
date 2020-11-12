clc; clear; close all;
global L R
L = 18e-3;
R = 2030;
[T,Res]=lyapunov(3,@ChuasCircuitSimFunc2,@ode45,0,0.00001,0.01,[0.01 -0.01 0],100);
plot(T,Res);
title('Dynamics of Lyapunov exponents');
xlabel('Time'); ylabel('Lyapunov exponents');

