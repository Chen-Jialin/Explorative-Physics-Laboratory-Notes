clc; clear; close all;
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
global R
R = 1800;

global Vmax1 Vmin1 Vmax2 Vmin2
Vmax1 = R3 / (R2 + R3) * V1;
Vmin1 = R3 / (R2 + R3) * V2;
Vmax2 = R6 / (R5 + R6) * V1;
Vmin2 = R6 / (R5 + R6) * V2;

% [T,Res]=lyapunov(3,@lorenz_ext,@ode45,0,5e-6,1e-2,[0.01 0.01 0],1000);
% plot(T,Res);
% xlim([0.002,0.01])
% grid on
% title('Dynamics of Lyapunov exponents');
% xlabel('Time'); ylabel('Lyapunov exponents');

R = 1950;
i = 0;
for V_C1 = -2:0.1:0
    i = i + 1;
    j = 0;
    for V_C2 = -2:0.1:2
        disp([V_C1, V_C2])
        j = j + 1;
        [T,Res]=lyapunov(3,@lorenz_ext,@ode45,0,5e-6,1e-2,[V_C1 V_C2 0],1000);
        LyapunovExponent(i, j) = Res(end, 1);
   end
end
LyapunovExponent(end + 1: 2 * end - 1, :) = LyapunovExponent(end - 1:-1:1, end:-1:1);
V_C1 = -2:0.1:2;
V_C2 = -2:0.1:2;
figure(1)
imagesc(V_C1, V_C2, log10(LyapunovExponent))
xlabel('V_{C1}(0) / V', 'fontsize', 16)
ylabel('V_{C2}(0) / V', 'fontsize', 16)
h = colorbar;
ylabel(h, 'log_{10}\lambda_{V_{C1}}', 'fontsize', 16)

LyapunovExponent_VC1 = zeros(size(1960:0.1:2000));
LyapunovExponent_VC2 = zeros(size(1960:0.1:2000));
LyapunovExponent_IL = zeros(size(1960:0.1:2000));
i = 1;
for R = 1960:0.1:2000
    disp(R)
    [T,Res]=lyapunov(3,@lorenz_ext,@ode45,0,5e-6,1e-2,[0.01 -0.01 0],1000);
    LyapunovExponent_VC1(i) = Res(end, 1);
    LyapunovExponent_VC2(i) = Res(end, 2);
    LyapunovExponent_IL(i) = Res(end, 3);
    i = i + 1;
end
figure(2)
subplot(3, 1, 1)
plot(1960:0.1:2000, LyapunovExponent_VC1)
grid on
set(gca,'xticklabel',{[]})
xline(1543.5, '--')
xline(1934.5, '--')
xline(2020.5, '--')
ylabel('V_{C1}的最大李雅普诺夫指数')
subplot(3, 1, 2)
plot(1960:0.1:2000, LyapunovExponent_VC2)
grid on
set(gca,'xticklabel',{[]})
xline(1543.5, '--')
xline(1934.5, '--')
xline(2020.5, '--')
ylabel('V_{C2}的最大李雅普诺夫指数')
subplot(3, 1, 3)
plot(1960:0.1:2000, LyapunovExponent_IL)
xlabel('R / \Omega', 'fontsize', 16)
ylabel('I_L的最大李雅普诺夫指数')
grid on
xline(1543.5, '--')
xline(1934.5, '--')
xline(2020.5, '--')