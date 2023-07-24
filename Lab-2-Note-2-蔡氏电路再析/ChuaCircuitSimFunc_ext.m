function f=ChuaCircuitSimFunc_ext(t,X)
%
%  Lorenz equation 
%
%               dx/dt = SIGMA*(y - x)
%               dy/dt = R*x - y -x*z
%               dz/dt= x*y - BETA*z
%
%        In demo run SIGMA = 10, R = 28, BETA = 8/3
%        Initial conditions: x(0) = 0, y(0) = 1, z(0) = 0;
%        Reference values for t=10 000 : 
%              L_1 = 0.9022, L_2 = 0.0003, LE3 = -14.5691
%
%        See:
%    K. Ramasubramanian, M.S. Sriram, "A comparative study of computation 
%    of Lyapunov spectra with different algorithms", Physica D 139 (2000) 72-86.
%
% --------------------------------------------------------------------
% Copyright (C) 2004, Govorukhin V.N.


% Values of parameters
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
L = 18e-3;
global R;

Vmax1 = R3 / (R2 + R3) * V1;
Vmin1 = R3 / (R2 + R3) * V2;
Vmax2 = R6 / (R5 + R6) * V1;
Vmin2 = R6 / (R5 + R6) * V2;

x=X(1); y=X(2); z=X(3);

I_R = (x <= Vmax1) .* (x >= Vmin1) .* (- R2 / (R3 * R1)) * x...
    + (x > Vmax1) .* (x - V1) / R1...
    + (x < Vmin1) .* (x - V2) / R1...
    + (x <= Vmax2) .* (x >= Vmin2) .* (- R5 / (R6 * R4)) * x...
    + (x > Vmax2) .* (x - V1) / R4...
    + (x < Vmin2) .* (x - V2) / R4;

dI_R = (x <= Vmax1) .* (x >= Vmin1) * (- R2 / (R3 * R1))...
    + (x > Vmax1) / R1...
    + (x < Vmin1) / R1...
    + (x <= Vmax2) .* (x >= Vmin2) * (- R5 / (R6 * R4))...
    + (x > Vmax2) / R4...
    + (x < Vmin2) / R4;

Y= [X(4), X(7), X(10);
    X(5), X(8), X(11);
    X(6), X(9), X(12)];

f=zeros(9,1);

%Lorenz equation
f(1) = - (I_R - (y - x) / R) / C1;
f(2) = - ((y - x) / R - z) / C2;
f(3) = - y / L;

%Linearized system

 Jac=[- dI_R / C1 - 1 / R / C1, 1 / R / C1, 0;
     1 / R / C2, - 1 / R / C2, 1 ./ C2;
     0, - 1 / L, 0];
  
%Variational equation   
f(4:12)=Jac*Y;

%Output data must be a column vector


