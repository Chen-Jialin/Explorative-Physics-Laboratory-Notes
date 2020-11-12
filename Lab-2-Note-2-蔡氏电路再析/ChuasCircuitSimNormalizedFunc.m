function dy = ChuasCircuitSimNormalizedFunc(t, y)
global beta

alpha = 10;
m0 = -0.1364;
m1 = 0.3865;

h = m1 * y(1, 1) + 0.5 * (m0 - m1) * (abs(y(1, 1) + 1) - abs(y(1, 1) - 1));

dy(1, 1) = alpha * (y(2, 1) - h);
dy(2, 1) = y(1, 1) - y(2, 1) + y(3, 1);
dy(3, 1) = - beta * y(2, 1);