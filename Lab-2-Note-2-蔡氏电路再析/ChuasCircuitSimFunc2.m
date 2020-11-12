function f = ChuasCircuitSimFunc2(t, y)
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

Y = [y(4), y(7), y(10);
    y(5), y(8), y(11);
    y(6), y(9), y(12)];

I_R = @(V) (V <= Vmax1) .* (V >= Vmin1) .* (- R2 / (R3 * R1) * V)...
    + (V > Vmax1) .* (V - V1) / R1...
    + (V < Vmin1) .* (V - V2) / R1...
    + (V <= Vmax2) .* (V >= Vmin2) .* (- R5 / (R6 * R4) * V)...
    + (V > Vmax2) .* (V - V1) / R4...
    + (V < Vmin2) .* (V - V2) / R4;
f(1,1) = - (I_R(y(1, 1)) - (y(2, 1) - y(1, 1)) / R) / C1;
f(2,1) = - ((y(2, 1) - y(1, 1)) / R - y(3, 1)) / C2;
f(3,1) = - y(2, 1) / L;

Jac=[-(y(1, 1) <= Vmax1) .* (y(1, 1) >= Vmin1) .* (- R2 / (R3 * R1))...
    - (y(1, 1) > Vmax1) .* 1 / R1...
    - (y(1, 1) < Vmin1) .* 1 / R1...
    - (y(1, 1) <= Vmax2) .* (y(1, 1) >= Vmin2) .* (- R5 / (R6 * R4))...
    - (y(1, 1) > Vmax2) .* 1 / R4...
    - (y(1, 1) < Vmin2) .* 1 / R4 - 1 / R / C1, 1 / R / C1, 0;
    1 / R / C2, - 1 / R / C2, 1 / C2;
    0, - 1 / L, 0];

%Variational equation   
f(4:12)=Jac * Y;