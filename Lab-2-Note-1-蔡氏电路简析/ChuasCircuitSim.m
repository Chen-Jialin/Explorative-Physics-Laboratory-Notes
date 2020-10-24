clc; clear; close all;
R1 = 220;
R2 = 220;
R3 = 2200;
R4 = 22000;
R5 = 22000;
R6 = 3300;
V1 = 9;
V2 = 9;
L1 = 18e-3;
C1 = 10e-9;
C2 = 100e-9;
global R; % Variable resistance
R = 1540;

% Negative impedance converter 1 I-V
Vmax1 = R3 / (R2 + R3) * V2;
Vmin1 = R3 / (R2 + R3) * (-V1);
V = linspace(2 * Vmin1, 2 * Vmax1, 101);
I = (V <= Vmax1) .* (V >= Vmin1) .* (- R2 / (R3 * R1) * V)...
    + (V > Vmax1) .* (V - V1) / R1...
    + (V < Vmin1) .* (V + V2) / R1;
figure(1)
plot(V, 1e3 * I, 'k-')
hold on

% Negative impedance converter 2 I-V
Vmax2 = R6 / (R5 + R6) * V2;
Vmin2 = R6 / (R5 + R6) * (-V1);
I = (V <= Vmax2) .* (V >= Vmin2) .* (- R5 / (R6 * R4) * V)...
    + (V > Vmax2) .* (V - V1) / R4...
    + (V < Vmin2) .* (V + V2) / R4;
plot(V, 1e3 * I, 'k--')

grid on
xlabel('V / V', 'fontsize', 16)
ylabel('I / mA', 'fontsize', 16)
legend('蔡氏二极管1（含R1,R2,R3）', '蔡氏二极管2（含R4,R5,R6）', 'fontsize', 12)

% Chua's Diode I-V
I = (V <= Vmax1) .* (V >= Vmin1) .* (- R2 / (R3 * R1) * V)...
    + (V > Vmax1) .* (V - V1) / R1...
    + (V < Vmin1) .* (V + V2) / R1...
    + (V <= Vmax2) .* (V >= Vmin2) .* (- R5 / (R6 * R4) * V)...
    + (V > Vmax2) .* (V - V1) / R4...
    + (V < Vmin2) .* (V + V2) / R4;
figure(2)
plot(V, 1e3 * I, 'k', 'linewidth', 2)
grid on
xlabel('V / V', 'fontsize', 16)
ylabel('I / mA', 'fontsize', 16)

% Chua's circuit simulation
% Simulation time range
t_span = [0, 0.05];
% Initial condition
y_0 = [0; 0.01; -0.02];
% y = [I_L; V_C1; V_C2]^T

% R = 1538;
% y_0 = [0; 0.01; -0.01];
% for i = 0:5
%     [t, y] = ode45(@ChuaCircuitSimFunc, t_span, y_0);
%     figure(3 + i)
%     plot(y(100:end,2), y(100:end,3), '-')
%     R = R + 1;
% end

R = 1540;
y = [0; 0.01; -0.02];
for i = 0:5
    [t, y] = ode45(@ChuaCircuitSimFunc, t_span, y_0);
    figure(3 + i)
    plot(y(:,2), y(:,3), '-')
    y_0(3) = y_0(3) + 0.01;
end