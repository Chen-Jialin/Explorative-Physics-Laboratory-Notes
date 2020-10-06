clc; clear; close all;
kDeltaL = -4 * pi: 0.01: 4 * pi;
I_I0 = (1 + cos(kDeltaL)) / 2;
plot(kDeltaL, I_I0);
xlim([-4 * pi, 4 * pi]);
xlabel('k\Delta{}L / rad', 'fontsize', 16);
ylabel('I / I_0', 'fontsize', 16);
