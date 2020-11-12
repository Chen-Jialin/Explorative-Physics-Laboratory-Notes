clc; clear; close all;
global beta
beta = 19;
for beta = 16:0.01:20
    disp(beta)
    t_span = [0, 500];
    y_0 = [0.2; 0.2; 0.2];
    [t, y] = ode45(@ChuasCircuitSimNormalizedFunc, t_span, y_0);
    y = y(t > 400, :);
    t = t(t > 400, 1);
    tmp = y(round(end / 1) - 500:end - 1, 1);
    tmp = tmp(tmp > y(round(end / 1) - 500 - 1: end - 2,1) & tmp > y(round(end / 1) - 500 + 1:end, 1));% & tmp > fixedpoint(2, 1));
    plot(beta, tmp, 'k.')
    hold on
end
xlabel('\beta', 'fontsize', 16)
ylabel('x', 'fontsize', 16)
% plot(y(:, 1), y(:, 2))