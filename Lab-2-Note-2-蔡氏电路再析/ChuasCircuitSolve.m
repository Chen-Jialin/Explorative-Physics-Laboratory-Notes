% 请选择需要的节运行
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
global L R 
L = 18e-3;
R = 1800;
global Vmax1 Vmin1 Vmax2 Vmin2
Vmax1 = R3 / (R2 + R3) * V1;
Vmin1 = R3 / (R2 + R3) * V2;
Vmax2 = R6 / (R5 + R6) * V1;
Vmin2 = R6 / (R5 + R6) * V2;
global turningpoint
global A1 b1 A2 b2 A3 b3 A4 b4 A5 b5

%% test
t_span = [0, 0.1];
y_0 = [0.01; -0.01; 0];
[t, y] = ode45(@ChuasCircuitSimFunc, t_span, y_0);
plot(y(round(end / 2):round(end / 1), 1), y(round(end / 2):round(end / 1), 2))

%% V_{C1}-V_{C2} evolution trace under different R
figure(1)
R = 10;
fixedpoint = fixedPoint();
t_span = [0, 0.01];
y_0 = [0.01; -0.01; 0];
% [t, y] = ode45(@ChuasCircuitSimFunc, t_span, y_0);
% plot(y(round(end / 2):round(end / 1), 1), y(round(end / 2):round(end / 1), 2))
% hold on
% % plot(fixedpoint(:, 1), fixedpoint(:, 2), '.' ,'markersize', 8)
% % xline(turningpoint(1), '-', 'linewidth', 2);
% % xline(turningpoint(2), '-', 'linewidth', 2);
% % xline(turningpoint(3), '-', 'linewidth', 2);
% % xline(turningpoint(4), '-', 'linewidth', 2);
% hold off
% set(gca, 'fontSize', 14);
% xlabel('V_{C1} / V', 'fontsize', 16)
% ylabel('V_{C2} / V', 'fontsize', 16)
% title(['R = ', num2str(R), ' \Omega'], 'fontsize', 16)
% grid on
% axis equal
% xlim([-11, 11])
% ylim([-11, 11])
for R = [100 500 1000 1500 1600 1700 1800 1900 1950 1965 1967.5 1970 1980 1990 1994 1998 2004 2006 2008]
    fixedpoint = fixedPoint();
    t_span = [0, 0.2];
    y_0 = [0.01; -0.01; 0];
    [t, y] = ode45(@ChuasCircuitSimFunc, t_span, y_0);
    plot(y(round(end / 2):round(end / 1), 1), y(round(end / 2):round(end / 1), 2))
    hold off
    set(gca, 'fontSize', 14);
    xlabel('V_{C1} / V', 'fontsize', 16)
    ylabel('V_{C2} / V', 'fontsize', 16)
    title(['R = ', num2str(R), ' \Omega'], 'fontsize', 16)
    grid on
end
R = 2100;
fixedpoint = fixedPoint();
t_span = [0, 0.2];
y_0 = [0.01; -0.01; 0];
[t, y] = ode45(@ChuasCircuitSimFunc, t_span, y_0);
plot(y(100:round(end / 1), 1), y(100:round(end / 1), 2))
hold off
xlabel('V_{C1} / V', 'fontsize', 16)
ylabel('V_{C2} / V', 'fontsize', 16)
title(['R = ', num2str(R), ' \Omega'], 'fontsize', 16)
grid on

% pic_num = 1;
% for R = 0:100:2000
%     [t, y] = ode45(@ChuasCircuitSimFunc, t_span, y_0);
%     plot(y(1:end, 1), y(1:end, 2))
%     title(['R = ', num2str(R)])
%     grid on
%     drawnow
%     F=getframe(gcf);
%     I=frame2im(F);
%     [I,map]=rgb2ind(I,256);
% 
%     if pic_num == 1
%     imwrite(I,map,'evolution.gif','gif','Loopcount',inf,'DelayTime',0.01);
% 
%     else
%     imwrite(I,map,'evolution.gif','gif','WriteMode','append','DelayTime',0.01);
% 
%     end
% 
%     pic_num = pic_num + 1;
% 
% end

%% max[V_{C1}] evolution under different L
R = 1800;
for L = 14.5e-3:0.01e-3:15.5e-3
    disp(L)
    t_span = [0, 0.1];
    y_0 = [0.01; -0.01; 0];
    [t, y] = ode45(@ChuasCircuitSimFunc, t_span, y_0);
    tmp = y(round(end / 2):end - 1, 1);
    tmp = tmp(tmp > y(round(end / 2) - 1: end - 2,1) & tmp > y(round(end / 2) + 1:end, 1));% & tmp > fixedpoint(2, 1));
    plot(L, tmp, 'k.')
    hold on
end
hold off

%% max[V_{C1}] evolution under different R - ode45
L = 18e-3;
for R = 1960:0.1:2000
% for R = 1920:1:2020
    disp(R)
    t_span = [0, 0.4];
    y_0 = [0.01; -0.01; 0];
    [t, y] = ode45(@ChuasCircuitSimFunc, t_span, y_0);
    tmp = y(round(end / 1) - 1000:end - 1, 1);
    tmp = tmp(tmp > y(round(end / 1) - 1000 - 1: end - 2, 1) & tmp > y(round(end / 1) - 1000 + 1:end, 1));% & tmp > fixedpoint(2, 1));
    plot(R, tmp, 'k.')
    hold on
end
xlabel('R / \Omega', 'fontsize', 16)
ylabel('V_{C1}极大值 / V', 'fontsize', 16)
hold off

%% sensitivity of the evolution trace to initial condition
R = 1950;
t_span = [0, 0.2];
y_0 = [0.01; -0.01; 0];
[t, y] = ode45(@ChuasCircuitSimFunc, t_span, y_0);
plot(y(:, 1), y(:, 2))
set(gca, 'fontSize', 14);
xlabel('V_{C1} / V', 'fontsize', 20)
ylabel('V_{C2} / V', 'fontsize', 20)
% title(['R = ', num2str(R), ' \Omega'], 'fontsize', 24)
xlim([-7, 7])
grid on

%% three views of evolution trace & fixed points
R = 1500;
fixedpoint = fixedPoint();
t_span = [0, 0.2];
y_0 = [0.01; -0.01; 0];
[t, y] = ode45(@ChuasCircuitSimFunc, t_span, y_0);

figure(2)
subplot(2, 2, 1)
plot(y(:, 1), y(:, 2), '-')
hold on
xlabel('V_{C1} / V')
ylabel('V_{C2} / V')
plot(fixedpoint(1:3, 1), fixedpoint(1:3, 2), '.' ,'markersize', 8)
hold off
grid on

subplot(2, 2, 2)
plot(1e3 * y(:, 3), y(:, 2), '-')
hold on
xlabel('I_L / mA')
ylabel('V_{C2} / V')
plot(1e3 * fixedpoint(1:3, 3), fixedpoint(1:3, 2), '.' ,'markersize', 8)
hold off
grid on

subplot(2, 2, 3)
plot(y(:, 1), 1e3 * y(:, 3), '-')
hold on
xlabel('V_{C1} / V')
ylabel('I_L / mA')
plot(fixedpoint(1:3, 1), 1e3 * fixedpoint(1:3, 3), '.' ,'markersize', 8)
hold off
grid on

attraction = repmat(y(:, 1) >= turningpoint(1), 1, 3) .* repmat(y(:, 1) <= turningpoint(2), 1, 3) .* (y * (A1^2).' + b1.' * A1.')...
    + repmat(y(:, 1) > turningpoint(2), 1, 3) .* repmat(y(:, 1) <= turningpoint(4), 1, 3) .* (y * (A2^2).' + b2.' * A2.')...
    + repmat(y(:, 1) < turningpoint(1), 1, 3) .* repmat(y(:, 1) >= turningpoint(3), 1, 3) .* (y * (A3^2).' + b3.' * A3.')...
    + repmat(y(:, 1) > turningpoint(4), 1, 3) .* (y * (A4^2).' + b4.' * A4.')...
    + repmat(y(:, 1) < turningpoint(3), 1, 3) .* (y * (A5^2).' + b5.' * A5.');

figure(3)
subplot(2, 2, 1)
plot(y(round(end / 2):end, 1), y(round(end / 2):end, 2), '-')
hold on
xlabel('V_{C1} / V', 'fontsize', 16)
ylabel('V_{C2} / V', 'fontsize', 16)
plot(fixedpoint(1:5, 1), fixedpoint(1:5, 2), '.' ,'markersize', 8)
quiver(y(round(end / 2):round(end / 1000):end, 1), y(round(end / 2):round(end / 1000):end, 2),...
    attraction(round(end / 2):round(end / 1000):end, 1), attraction(round(end / 2):round(end / 1000):end, 2), 2, 'g-', 'linewidth', 1)
xline(turningpoint(1), '-', 'linewidth', 2);
xline(turningpoint(2), '-', 'linewidth', 2);
xline(turningpoint(3), '-', 'linewidth', 2);
xline(turningpoint(4), '-', 'linewidth', 2);
hold off
grid on

subplot(2, 2, 2)
plot(1e3 * y(round(end / 2):end, 3), y(round(end / 2):end, 2), '-')
hold on
xlabel('I_L / mA', 'fontsize', 16)
ylabel('V_{C2} / V', 'fontsize', 16)
plot(1e3 * fixedpoint(1:5, 3), fixedpoint(1:5, 2), '.' ,'markersize', 8)
quiver(1e3 * y(round(end / 2):round(end / 1000):end, 3), y(round(end / 2):round(end / 1000):end, 2),...
    1e3 * attraction(round(end / 2):round(end / 1000):end, 3), attraction(round(end / 2):round(end / 1000):end, 2), 2, 'g-', 'linewidth', 1)
hold off
grid on

subplot(2, 2, 3)
plot(y(round(end / 2):end, 1), 1e3 * y(round(end / 2):end, 3), '-')
hold on
xlabel('V_{C1} / V', 'fontsize', 16)
ylabel('I_L / mA', 'fontsize', 16)
plot(fixedpoint(1:5, 1), 1e3 * fixedpoint(1:5, 3), '.' ,'markersize', 8)
quiver(y(round(end / 2):round(end / 1000):end, 1), 1e3 * y(round(end / 2):round(end / 1000):end, 3),...
    attraction(round(end / 2):round(end / 1000):end, 1), 1e3 * attraction(round(end / 2):round(end / 1000):end, 3), 2, 'g-', 'linewidth', 1)
xline(turningpoint(1), '-', 'linewidth', 2);
xline(turningpoint(2), '-', 'linewidth', 2);
xline(turningpoint(3), '-', 'linewidth', 2);
xline(turningpoint(4), '-', 'linewidth', 2);
hold off
grid on

figure(4)
i = 1;
V_C1 = -5:0.25:5;
V_C2 = -1.5:0.5:1.5;
[V_C1, V_C2] = meshgrid(V_C1, V_C2);
V_C1 = reshape(V_C1, [], 1);
V_C2 = reshape(V_C2, [], 1);
for I_L = 1e-3 * (12:-3:-12)
    subplot(9, 1, i)
    attraction = (V_C1 >= turningpoint(1)) .* (V_C1 <= turningpoint(2)) .* ([V_C1, V_C2, I_L * ones(size(V_C1))] * (A1^2).' + b1.' * A1.')...
        + (V_C1 > turningpoint(2)) .* (V_C1 <= turningpoint(4)) .* ([V_C1, V_C2, I_L * ones(size(V_C1))] * (A2^2).' + b2.' * A2.')...
        + (V_C1 < turningpoint(1)) .* (V_C1 >= turningpoint(3)) .* ([V_C1, V_C2, I_L * ones(size(V_C1))] * (A3^2).' + b3.' * A3.')...
        + (V_C1 > turningpoint(4)) .* ([V_C1, V_C2, I_L * ones(size(V_C1))] * (A4^2).' + b4.' * A4.')...
        + (V_C1 < turningpoint(3)) .* ([V_C1, V_C2, I_L * ones(size(V_C1))] * (A5^2).' + b5.' * A5.');
    imagesc(-5:0.25:5, -1.5:0.5:1.5, reshape(attraction(:, 3), length(-5:0.25:5), []))
    colormap hot
    h = colorbar;
    caxis([-8e6, 8d6])
    hold on
    quiver(V_C1, V_C2, attraction(:, 1), attraction(:, 2), 'g-', 'linewidth', 1)
    xline(turningpoint(1), '-', 'linewidth', 2);
    xline(turningpoint(2), '-', 'linewidth', 2);
    hold off
    if (I_L == -12e-3)
        xlabel('V_{C1} / V', 'fontsize', 16)
    end
    if (I_L == 0)
        ylabel('V_{C2} / V', 'fontsize', 16)
        ylabel(h, '吸引力在I_L轴方向的分量', 'fontsize', 16)
    end
    i = i + 1;
end

%% three views of evolution trace & fixed points
R = 1800;
fixedpoint = fixedPoint();
t_span = [0, 0.2];
y_0 = [0.01; -0.01; 0];
[t, y] = ode45(@ChuasCircuitSimFunc, t_span, y_0);

figure(2)
subplot(2, 2, 1)
plot(y(:, 1), y(:, 2), '-')
hold on
xlabel('V_{C1} / V')
ylabel('V_{C2} / V')
plot(fixedpoint(1:3, 1), fixedpoint(1:3, 2), '.' ,'markersize', 8)
hold off
grid on

subplot(2, 2, 2)
plot(1e3 * y(:, 3), y(:, 2), '-')
hold on
xlabel('I_L / mA')
ylabel('V_{C2} / V')
plot(1e3 * fixedpoint(1:3, 3), fixedpoint(1:3, 2), '.' ,'markersize', 8)
hold off
grid on

subplot(2, 2, 3)
plot(y(:, 1), 1e3 * y(:, 3), '-')
hold on
xlabel('V_{C1} / V')
ylabel('I_L / mA')
plot(fixedpoint(1:3, 1), 1e3 * fixedpoint(1:3, 3), '.' ,'markersize', 8)
hold off
grid on

attraction = repmat(y(:, 1) >= turningpoint(1), 1, 3) .* repmat(y(:, 1) <= turningpoint(2), 1, 3) .* (y * (A1^2).' + b1.' * A1.')...
    + repmat(y(:, 1) > turningpoint(2), 1, 3) .* repmat(y(:, 1) <= turningpoint(4), 1, 3) .* (y * (A2^2).' + b2.' * A2.')...
    + repmat(y(:, 1) < turningpoint(1), 1, 3) .* repmat(y(:, 1) >= turningpoint(3), 1, 3) .* (y * (A3^2).' + b3.' * A3.')...
    + repmat(y(:, 1) > turningpoint(4), 1, 3) .* (y * (A4^2).' + b4.' * A4.')...
    + repmat(y(:, 1) < turningpoint(3), 1, 3) .* (y * (A5^2).' + b5.' * A5.');

figure(3)
subplot(2, 2, 1)
plot(y(round(end / 2):end, 1), y(round(end / 2):end, 2), '-')
hold on
xlabel('V_{C1} / V', 'fontsize', 16)
ylabel('V_{C2} / V', 'fontsize', 16)
plot(fixedpoint(1:5, 1), fixedpoint(1:5, 2), '.' ,'markersize', 8)
quiver(y(round(end / 2):round(end / 1000):end, 1), y(round(end / 2):round(end / 1000):end, 2),...
    attraction(round(end / 2):round(end / 1000):end, 1), attraction(round(end / 2):round(end / 1000):end, 2), 2, 'g-', 'linewidth', 1)
xline(turningpoint(1), '-', 'linewidth', 2);
xline(turningpoint(2), '-', 'linewidth', 2);
xline(turningpoint(3), '-', 'linewidth', 2);
xline(turningpoint(4), '-', 'linewidth', 2);
hold off
grid on

subplot(2, 2, 2)
plot(1e3 * y(round(end / 2):end, 3), y(round(end / 2):end, 2), '-')
hold on
xlabel('I_L / mA', 'fontsize', 16)
ylabel('V_{C2} / V', 'fontsize', 16)
plot(1e3 * fixedpoint(1:5, 3), fixedpoint(1:5, 2), '.' ,'markersize', 8)
quiver(1e3 * y(round(end / 2):round(end / 1000):end, 3), y(round(end / 2):round(end / 1000):end, 2),...
    1e3 * attraction(round(end / 2):round(end / 1000):end, 3), attraction(round(end / 2):round(end / 1000):end, 2), 2, 'g-', 'linewidth', 1)
hold off
grid on

subplot(2, 2, 3)
plot(y(round(end / 2):end, 1), 1e3 * y(round(end / 2):end, 3), '-')
hold on
xlabel('V_{C1} / V', 'fontsize', 16)
ylabel('I_L / mA', 'fontsize', 16)
plot(fixedpoint(1:5, 1), 1e3 * fixedpoint(1:5, 3), '.' ,'markersize', 8)
quiver(y(round(end / 2):round(end / 1000):end, 1), 1e3 * y(round(end / 2):round(end / 1000):end, 3),...
    attraction(round(end / 2):round(end / 1000):end, 1), 1e3 * attraction(round(end / 2):round(end / 1000):end, 3), 2, 'g-', 'linewidth', 1)
xline(turningpoint(1), '-', 'linewidth', 2);
xline(turningpoint(2), '-', 'linewidth', 2);
xline(turningpoint(3), '-', 'linewidth', 2);
xline(turningpoint(4), '-', 'linewidth', 2);
hold off
grid on

figure(4)
i = 1;
V_C1 = -5:0.25:5;
V_C2 = -1.5:0.5:1.5;
[V_C1, V_C2] = meshgrid(V_C1, V_C2);
V_C1 = reshape(V_C1, [], 1);
V_C2 = reshape(V_C2, [], 1);
for I_L = 1e-3 * (4:-1:-4)
    subplot(9, 1, i)
    attraction = (V_C1 >= turningpoint(1)) .* (V_C1 <= turningpoint(2)) .* ([V_C1, V_C2, I_L * ones(size(V_C1))] * (A1^2).' + b1.' * A1.')...
        + (V_C1 > turningpoint(2)) .* (V_C1 <= turningpoint(4)) .* ([V_C1, V_C2, I_L * ones(size(V_C1))] * (A2^2).' + b2.' * A2.')...
        + (V_C1 < turningpoint(1)) .* (V_C1 >= turningpoint(3)) .* ([V_C1, V_C2, I_L * ones(size(V_C1))] * (A3^2).' + b3.' * A3.')...
        + (V_C1 > turningpoint(4)) .* ([V_C1, V_C2, I_L * ones(size(V_C1))] * (A4^2).' + b4.' * A4.')...
        + (V_C1 < turningpoint(3)) .* ([V_C1, V_C2, I_L * ones(size(V_C1))] * (A5^2).' + b5.' * A5.');
    imagesc(-5:0.25:5, -1.5:0.5:1.5, reshape(attraction(:, 3), length(-5:0.25:5), []))
    colormap hot
    h = colorbar;
    caxis([-4e6, 4d6])
    hold on
    quiver(V_C1, V_C2, attraction(:, 1), attraction(:, 2), 'g-', 'linewidth', 1)
    xline(turningpoint(1), '-', 'linewidth', 2);
    xline(turningpoint(2), '-', 'linewidth', 2);
    hold off
    if (I_L == -4e-3)
        xlabel('V_{C1} / V', 'fontsize', 16)
    end
    if (I_L == 0)
        ylabel('V_{C2} / V', 'fontsize', 16)
        ylabel(h, '吸引力在I_L轴方向的分量', 'fontsize', 16)
    end
    i = i + 1;
end

%% three views of evolution trace & fixed points
R = 1950;
fixedpoint = fixedPoint();
t_span = [0, 0.2];
y_0 = [0.01; -0.01; 0];
[t, y] = ode45(@ChuasCircuitSimFunc, t_span, y_0);

figure(2)
subplot(2, 2, 1)
plot(y(:, 1), y(:, 2), '-')
hold on
xlabel('V_{C1} / V', 'fontsize', 16)
ylabel('V_{C2} / V', 'fontsize', 16)
plot(fixedpoint(1:3, 1), fixedpoint(1:3, 2), '.' ,'markersize', 8)
hold off
grid on

subplot(2, 2, 2)
plot(1e3 * y(:, 3), y(:, 2), '-')
hold on
xlabel('I_L / mA', 'fontsize', 16)
ylabel('V_{C2} / V', 'fontsize', 16)
plot(1e3 * fixedpoint(1:3, 3), fixedpoint(1:3, 2), '.' ,'markersize', 8)
hold off
grid on

subplot(2, 2, 3)
plot(y(:, 1), 1e3 * y(:, 3), '-')
hold on
xlabel('V_{C1} / V', 'fontsize', 16)
ylabel('I_L / mA', 'fontsize', 16)
plot(fixedpoint(1:3, 1), 1e3 * fixedpoint(1:3, 3), '.' ,'markersize', 8)
hold off
grid on

attraction = repmat(y(:, 1) >= turningpoint(1), 1, 3) .* repmat(y(:, 1) <= turningpoint(2), 1, 3) .* (y * (A1^2).' + b1.' * A1.')...
    + repmat(y(:, 1) > turningpoint(2), 1, 3) .* repmat(y(:, 1) <= turningpoint(4), 1, 3) .* (y * (A2^2).' + b2.' * A2.')...
    + repmat(y(:, 1) < turningpoint(1), 1, 3) .* repmat(y(:, 1) >= turningpoint(3), 1, 3) .* (y * (A3^2).' + b3.' * A3.')...
    + repmat(y(:, 1) > turningpoint(4), 1, 3) .* (y * (A4^2).' + b4.' * A4.')...
    + repmat(y(:, 1) < turningpoint(3), 1, 3) .* (y * (A5^2).' + b5.' * A5.');

figure(3)
subplot(2, 2, 1)
plot(y(round(end / 2):end, 1), y(round(end / 2):end, 2), '-')
hold on
xlabel('V_{C1} / V', 'fontsize', 16)
ylabel('V_{C2} / V', 'fontsize', 16)
plot(fixedpoint(1:3, 1), fixedpoint(1:3, 2), '.' ,'markersize', 8)
quiver(y(round(end / 2):round(end / 1000):end, 1), y(round(end / 2):round(end / 1000):end, 2),...
    attraction(round(end / 2):round(end / 1000):end, 1), attraction(round(end / 2):round(end / 1000):end, 2), 2, 'g-', 'linewidth', 1)
xline(turningpoint(1), '-', 'linewidth', 2);
xline(turningpoint(2), '-', 'linewidth', 2);
hold off
grid on

subplot(2, 2, 2)
plot(1e3 * y(round(end / 2):end, 3), y(round(end / 2):end, 2), '-')
hold on
xlabel('I_L / mA', 'fontsize', 16)
ylabel('V_{C2} / V', 'fontsize', 16)
plot(1e3 * fixedpoint(1:3, 3), fixedpoint(1:3, 2), '.' ,'markersize', 8)
quiver(1e3 * y(round(end / 2):round(end / 1000):end, 3), y(round(end / 2):round(end / 1000):end, 2),...
    1e3 * attraction(round(end / 2):round(end / 1000):end, 3), attraction(round(end / 2):round(end / 1000):end, 2), 2, 'g-', 'linewidth', 1)
hold off
grid on

subplot(2, 2, 3)
plot(y(round(end / 2):end, 1), 1e3 * y(round(end / 2):end, 3), '-')
hold on
xlabel('V_{C1} / V', 'fontsize', 16)
ylabel('I_L / mA', 'fontsize', 16)
plot(fixedpoint(1:3, 1), 1e3 * fixedpoint(1:3, 3), '.' ,'markersize', 8)
quiver(y(round(end / 2):round(end / 1000):end, 1), 1e3 * y(round(end / 2):round(end / 1000):end, 3),...
    attraction(round(end / 2):round(end / 1000):end, 1), 1e3 * attraction(round(end / 2):round(end / 1000):end, 3), 2, 'g-', 'linewidth', 1)
xline(turningpoint(1), '-', 'linewidth', 2);
xline(turningpoint(2), '-', 'linewidth', 2);
hold off
grid on

figure(4)
i = 1;
V_C1 = -2:0.25:8;
V_C2 = -1.5:0.5:1.5;
[V_C1, V_C2] = meshgrid(V_C1, V_C2);
V_C1 = reshape(V_C1, [], 1);
V_C2 = reshape(V_C2, [], 1);
for I_L = 1e-3 * (4:-1:-4)
    subplot(9, 1, i)
    attraction = (V_C1 >= turningpoint(1)) .* (V_C1 <= turningpoint(2)) .* ([V_C1, V_C2, I_L * ones(size(V_C1))] * (A1^2).' + b1.' * A1.')...
        + (V_C1 > turningpoint(2)) .* (V_C1 <= turningpoint(4)) .* ([V_C1, V_C2, I_L * ones(size(V_C1))] * (A2^2).' + b2.' * A2.')...
        + (V_C1 < turningpoint(1)) .* (V_C1 >= turningpoint(3)) .* ([V_C1, V_C2, I_L * ones(size(V_C1))] * (A3^2).' + b3.' * A3.')...
        + (V_C1 > turningpoint(4)) .* ([V_C1, V_C2, I_L * ones(size(V_C1))] * (A4^2).' + b4.' * A4.')...
        + (V_C1 < turningpoint(3)) .* ([V_C1, V_C2, I_L * ones(size(V_C1))] * (A5^2).' + b5.' * A5.');
    imagesc(-2:0.25:8, -1.5:0.5:1.5, reshape(attraction(:, 3), length(-2:0.25:8), []))
    colormap hot
    h = colorbar;
    caxis([-4e6, 3d6])
    hold on
    quiver(V_C1, V_C2, attraction(:, 1), attraction(:, 2), 'g-', 'linewidth', 1)
    xline(turningpoint(1), '-', 'linewidth', 2);
    xline(turningpoint(2), '-', 'linewidth', 2);
    hold off
    if (I_L == -4e-3)
        xlabel('V_{C1} / V', 'fontsize', 16)
    end
    if (I_L == 0)
        ylabel('V_{C2} / V', 'fontsize', 16)
        ylabel(h, '吸引力在I_L轴方向的分量', 'fontsize', 16)
    end
    i = i + 1;
end
%% three views of evolution trace & fixed points
R = 2100;
fixedpoint = fixedPoint();
t_span = [0, 0.2];
y_0 = [0.01; -0.01; 0];
[t, y] = ode45(@ChuasCircuitSimFunc, t_span, y_0);

figure(2)
subplot(2, 2, 1)
plot(y(:, 1), y(:, 2), '-')
hold on
xlabel('V_{C1} / V', 'fontsize', 16)
ylabel('V_{C2} / V', 'fontsize', 16)
plot(fixedpoint(1:3, 1), fixedpoint(1:3, 2), '.' ,'markersize', 8)
hold off
grid on

subplot(2, 2, 2)
plot(1e3 * y(:, 3), y(:, 2), '-')
hold on
xlabel('I_L / mA', 'fontsize', 16)
ylabel('V_{C2} / V', 'fontsize', 16)
plot(1e3 * fixedpoint(1:3, 3), fixedpoint(1:3, 2), '.' ,'markersize', 8)
hold off
grid on

subplot(2, 2, 3)
plot(y(:, 1), 1e3 * y(:, 3), '-')
hold on
xlabel('V_{C1} / V', 'fontsize', 16)
ylabel('I_L / mA', 'fontsize', 16)
plot(fixedpoint(1:3, 1), 1e3 * fixedpoint(1:3, 3), '.' ,'markersize', 8)
hold off
grid on

attraction = repmat(y(:, 1) >= turningpoint(1), 1, 3) .* repmat(y(:, 1) <= turningpoint(2), 1, 3) .* (y * (A1^2).' + b1.' * A1.')...
    + repmat(y(:, 1) > turningpoint(2), 1, 3) .* repmat(y(:, 1) <= turningpoint(4), 1, 3) .* (y * (A2^2).' + b2.' * A2.')...
    + repmat(y(:, 1) < turningpoint(1), 1, 3) .* repmat(y(:, 1) >= turningpoint(3), 1, 3) .* (y * (A3^2).' + b3.' * A3.')...
    + repmat(y(:, 1) > turningpoint(4), 1, 3) .* (y * (A4^2).' + b4.' * A4.')...
    + repmat(y(:, 1) < turningpoint(3), 1, 3) .* (y * (A5^2).' + b5.' * A5.');

figure(3)
subplot(2, 2, 1)
plot(y(1:end, 1), y(1:end, 2), '-')
hold on
xlabel('V_{C1} / V', 'fontsize', 16)
ylabel('V_{C2} / V', 'fontsize', 16)
plot(fixedpoint(1:3, 1), fixedpoint(1:3, 2), '.' ,'markersize', 8)
quiver(y(1:round(end / 1000):end, 1), y(1:round(end / 1000):end, 2),...
    attraction(1:round(end / 1000):end, 1), attraction(1:round(end / 1000):end, 2), 2, 'g-', 'linewidth', 1)
xline(turningpoint(1), '-', 'linewidth', 2);
xline(turningpoint(2), '-', 'linewidth', 2);
hold off
grid on

subplot(2, 2, 2)
plot(1e3 * y(1:end, 3), y(1:end, 2), '-')
hold on
xlabel('I_L / mA', 'fontsize', 16)
ylabel('V_{C2} / V', 'fontsize', 16)
plot(1e3 * fixedpoint(1:3, 3), fixedpoint(1:3, 2), '.' ,'markersize', 8)
quiver(1e3 * y(1:round(end / 1000):end, 3), y(1:round(end / 1000):end, 2),...
    1e3 * attraction(1:round(end / 1000):end, 3), attraction(1:round(end / 1000):end, 2), 2, 'g-', 'linewidth', 1)
hold off
grid on

subplot(2, 2, 3)
plot(y(1:end, 1), 1e3 * y(1:end, 3), '-')
hold on
xlabel('V_{C1} / V', 'fontsize', 16)
ylabel('I_L / mA', 'fontsize', 16)
plot(fixedpoint(1:3, 1), 1e3 * fixedpoint(1:3, 3), '.' ,'markersize', 8)
quiver(y(1:round(end / 1000):end, 1), 1e3 * y(1:round(end / 1000):end, 3),...
    attraction(1:round(end / 1000):end, 1), 1e3 * attraction(1:round(end / 1000):end, 3), 2, 'g-', 'linewidth', 1)
xline(turningpoint(1), '-', 'linewidth', 2);
xline(turningpoint(2), '-', 'linewidth', 2);
hold off
grid on

figure(4)
i = 1;
V_C1 = -2:0.2:7;
V_C2 = -0.5:0.1:0.5;
[V_C1, V_C2] = meshgrid(V_C1, V_C2);
V_C1 = reshape(V_C1, [], 1);
V_C2 = reshape(V_C2, [], 1);
for I_L = 1e-3 * (4:-1:-4)
    subplot(9, 1, i)
    attraction = (V_C1 >= turningpoint(1)) .* (V_C1 <= turningpoint(2)) .* ([V_C1, V_C2, I_L * ones(size(V_C1))] * (A1^2).' + b1.' * A1.')...
        + (V_C1 > turningpoint(2)) .* (V_C1 <= turningpoint(4)) .* ([V_C1, V_C2, I_L * ones(size(V_C1))] * (A2^2).' + b2.' * A2.')...
        + (V_C1 < turningpoint(1)) .* (V_C1 >= turningpoint(3)) .* ([V_C1, V_C2, I_L * ones(size(V_C1))] * (A3^2).' + b3.' * A3.')...
        + (V_C1 > turningpoint(4)) .* ([V_C1, V_C2, I_L * ones(size(V_C1))] * (A4^2).' + b4.' * A4.')...
        + (V_C1 < turningpoint(3)) .* ([V_C1, V_C2, I_L * ones(size(V_C1))] * (A5^2).' + b5.' * A5.');
    imagesc(-2:0.2:7, -0.5:0.1:0.5, reshape(attraction(:, 3), length(-2:0.2:7), []))
    colormap hot
    h = colorbar;
    caxis([-4e6, 2.5d6])
    hold on
    quiver(V_C1, V_C2, attraction(:, 1), attraction(:, 2), 'g-', 'linewidth', 1)
    xline(turningpoint(1), '-', 'linewidth', 2);
    xline(turningpoint(2), '-', 'linewidth', 2);
    hold off
    if (I_L == -4e-3)
        xlabel('V_{C1} / V', 'fontsize', 16)
    end
    if (I_L == 0)
        ylabel('V_{C2} / V', 'fontsize', 16)
        ylabel(h, '吸引力在I_L轴方向的分量', 'fontsize', 16)
    end
    i = i + 1;
end
%% fixed points & turning points
function fixedpoints = fixedPoint()
% fixed points
R1 = 220;
R2 = 220;
R3 = 2200;
R4 = 22000;
R5 = 22000;
R6 = 3300;
V1 = 9;
V2 = -9;
L1 = 18e-3;
C1 = 10e-9;
C2 = 100e-9;
global R

% turning points
global turningpoint
turningpoint(1, 1) = R6 * V2 / (R5 + R6);
turningpoint(2, 1) = R6 * V1 / (R5 + R6);
turningpoint(3, 1) = R3 * V2 / (R2 + R3);
turningpoint(4, 1) = R3 * V1 / (R2 + R3);

% R6 * Vnn / (R5 + R6) <= VC1 <= R6 * Vpp / (R5 + R6)
global A1
A1 = [- 1 / C1 / R + R2 / R3 / R1 / C1 + R5 / R6 / R4 / C1, 1 / C1 / R, 0;...
    1 / C2 / R, - 1 / C2 / R, 1 / C2;...
    0, - 1 / L1, 0];
global b1
b1 = [0; 0; 0];
fixedpoints(1, :) = - (A1 \ b1).';

% R6 * Vpp / (R5 + R6) < VC1 <= R3 * Vpp / (R2 + R3)
global A2
A2 = [- 1 / C1 / R + R2 / R3 / R1 / C1 - 1 / R4 / C1, 1 / C1 / R, 0;...
    1 / C2 / R, -1 / C2 / R, 1 / C2;...
    0, -1 / L1, 0];
global b2
b2 = [V1 / C1 / R4; 0; 0];
fixedpoints(2, :) = - (A2 \ b2).';

% R3 * Vnn / (R2 + R3) <= VC1 < R6 * Vnn / (R5 + R6)
global A3
A3 = [- 1 / C1 / R + R2 / R3 / R1 / C1 - 1 / R4 / C1, 1 / C1 / R, 0;...
    1 / C2 / R, -1 / C2 / R, 1 / C2;...
    0, -1 / L1, 0];
global b3
b3 = [V2 / C1 / R4; 0; 0];
fixedpoints(3, :) = - (A3 \ b3).';

% VC1 > R3 * Vpp / (R2 + R3)
global A4
A4 = [- 1 / C1 / R - 1 / R1 / C1 - 1 / R4 / C1, 1 / C1 / R, 0;...
    1 / C2 / R, -1 / C2 / R, 1 / C2;...
    0, -1 / L1, 0];
global b4
b4 = [V1 / C1 / R1 + V1 / C1 / R4; 0; 0];
fixedpoints(4, :) = - (A4 \ b4).';

% VC1 < R3 * Vnn / (R2 + R3)
global A5
A5 = [- 1 / C1 / R - 1 / R1 / C1 - 1 / R4 / C1, 1 / C1 / R, 0;...
    1 / C2 / R, -1 / C2 / R, 1 / C2;...
    0, -1 / L1, 0];
global b5
b5 = [V2 / C1 / R1 + V2 / C1 / R4; 0; 0];
fixedpoints(5, :) = - (A5 \ b5).';
end