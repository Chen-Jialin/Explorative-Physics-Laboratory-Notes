clc;clear;close all;
Deltaomega_gamma = -20:0.01:20;

i = 1;
for tau0 = [0.4, 2, 20]
    figure(i);
    i = i + 1;
    hold on;
    for kDeltaL = (0:5) * pi / 5
        tau = tau0 ./ (Deltaomega_gamma.^2 + 1);
        delta = - Deltaomega_gamma .* tau;
        I_I0 = (1 + exp(-2 * tau) + 2 * exp(-tau) .* cos(kDeltaL + delta)) / 4;
        plot(Deltaomega_gamma, I_I0);
    end
    xlabel('\Delta\omega/\gamma', 'fontsize', 16);
    ylabel('I/I_0', 'fontsize', 16);
    legend('(k\Delta L) = 0', '(k\Delta L) = \pi/5', '(k\Delta L) = 2\pi/5', '(k\Delta L) = 3\pi/5', '(k\Delta L) = 4\pi/5', '(k\Delta L) = \pi');
end