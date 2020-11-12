function dy = ChuasCircuitSimFunc(t, y)
R1 = 220;
R2 = 220;
R3 = 2200;
R4 = 22000;
R5 = 22000;
R6 = 3300;
V1 = 9;
V2 = -9;
C1 = 10e-9;
C2 = 100e-9;
global R L
global Vmax1 Vmin1 Vmax2 Vmin2
% Vmax1 = R3 / (R2 + R3) * V1;
% Vmin1 = R3 / (R2 + R3) * V2;
% Vmax2 = R6 / (R5 + R6) * V1;
% Vmin2 = R6 / (R5 + R6) * V2;
I_R = (y(1, 1) <= Vmax1) .* (y(1, 1) >= Vmin1) .* (- R2 / (R3 * R1) * y(1, 1))...
    + (y(1, 1) > Vmax1) .* (y(1, 1) - V1) / R1...
    + (y(1, 1) < Vmin1) .* (y(1, 1) - V2) / R1...
    + (y(1, 1) <= Vmax2) .* (y(1, 1) >= Vmin2) .* (- R5 / (R6 * R4) * y(1, 1))...
    + (y(1, 1) > Vmax2) .* (y(1, 1) - V1) / R4...
    + (y(1, 1) < Vmin2) .* (y(1, 1) - V2) / R4;
dy(1,1) = - (I_R - (y(2, 1) - y(1, 1)) / R) / C1;
dy(2,1) = - ((y(2, 1) - y(1, 1)) / R - y(3, 1)) / C2;
dy(3,1) = - y(2, 1) / L;