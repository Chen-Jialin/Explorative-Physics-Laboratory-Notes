function dy = ChuaCircuitSimFunc(t, y)
    R1 = 220;
    R2 = 220;
    R3 = 2200;
    R4 = 22000;
    R5 = 22000;
    R6 = 3300;
    V1 = 9;
    V2 = 9;
    L = 18e-3;
    C1 = 10e-9;
    C2 = 100e-9;
    global R;
    Vmax1 = R3 / (R2 + R3) * V2;
    Vmin1 = R3 / (R2 + R3) * (-V1);
    Vmax2 = R6 / (R5 + R6) * V2;
    Vmin2 = R6 / (R5 + R6) * (-V1);
    I_R = @(V) (V <= Vmax1) .* (V >= Vmin1) .* (- R2 / (R3 * R1) * V)...
    + (V > Vmax1) .* (V - V1) / R1...
    + (V < Vmin1) .* (V + V2) / R1...
    + (V <= Vmax2) .* (V >= Vmin2) .* (- R5 / (R6 * R4) * V)...
    + (V > Vmax2) .* (V - V1) / R4...
    + (V < Vmin2) .* (V + V2) / R4;
    dy(1,1) = - y(3, 1) / L;
    dy(2,1) = - (I_R(y(2, 1)) - (y(3, 1) - y(2, 1)) / R) / C1;
    dy(3,1) = - ((y(3, 1) - y(2, 1)) / R - y(1, 1)) / C2;