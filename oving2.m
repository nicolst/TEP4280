%% Initial commands
clear all
close all
clc

%% Computation of trajectory for $\Delta t = 0.001$ (a)

% Time step, end time and time steps
step = 0.001;
endt = 10.0;
N = ceil(endt / step);
step = endt / N;

% Magnet "radius"
global R
R = 1.0;

% Pendulum "height"
global d
d = 0.1*R;

% Initial conditions
x0 = 2.0;
y0 = 0.3;
u0 = 0.0;
v0 = 0.0;

% Define trajectory arrays
euler = zeros(2, N);
heun = zeros(2, N);

euler(1, 1) = x0;
euler(2, 1) = y0;
eY = [x0; y0; u0; v0];

heun(1, 1) = x0;
heun(2, 1) = y0;
hY = [x0; y0; u0; v0];

% Compute trajectories for $\Delta t = 0.001$
for i=2:N
    % Euler
    k = f(0, eY);
    eY = eY + step * k;
    euler(:, i) = eY(1:2);

    % Heun
    k1 = f(0, hY);
    k2 = f(0, hY + step * k1);
    hY = hY + (step / 2) * (k1 + k2);
    heun(:, i) = hY(1:2);
end

% ode45
[ot, oy] = ode45(@f, [0 endt], [x0; y0; u0; v0]);

%% Plot for $\Delta t = 0.001$ (b)
% We see from the plot that ode45 and Heun's method seem to coincide for
% quite som time before diverging, while Euler's method produces a wildly
% different result.
plot(heun(1,:), heun(2,:), ...
     euler(1,:), euler(2,:), '--', ...
     oy(:,1), oy(:,2), '-.')
legend({'Heun', 'Euler', 'ode45'}, 'Location', 'best')
xlabel('x')
ylabel('y')
title('Trajectory for pendulum, \Deltat = 0.001')
snapnow

%% Compute and plot trajectories for $\Delta t = 0.0001$ (c)

% Set new time step and reset initial conditions
step = 0.0001;
N = ceil(endt / step);
step = endt / N;

euler(1, 1) = x0;
euler(2, 1) = y0;
eY = [x0; y0; u0; v0];

heun(1, 1) = x0;
heun(2, 1) = y0;
hY = [x0; y0; u0; v0];

% Compute trajectory
for i=2:N
    % Euler
    k = f(0, eY);
    eY = eY + step * k;
    euler(:, i) = eY(1:2);

    % Heun
    k1 = f(0, hY);
    k2 = f(0, hY + step * k1);
    hY = hY + (step / 2) * (k1 + k2);
    heun(:, i) = hY(1:2);
end

% Plot trajectory
plot(heun(1,:), heun(2,:), ...
     euler(1,:), euler(2,:), '--', ...
     oy(:,1), oy(:,2), '-.')
legend({'Heun', 'Euler', 'ode45'}, 'Location', 'best')
xlabel('x')
ylabel('y')
title('Trajectory for pendulum, \Deltat = 0.0001')
snapnow

%% Compute and plot trajectory for $\Delta t = 0.00001$ (c)

% Set new time step and reset initial conditions
step = 0.00001;
N = ceil(endt / step);
step = endt / N;

euler(1, 1) = x0;
euler(2, 1) = y0;
eY = [x0; y0; u0; v0];

heun(1, 1) = x0;
heun(2, 1) = y0;
hY = [x0; y0; u0; v0];

% Compute trajectory
for i=2:N
    % Euler
    k = f(0, eY);
    eY = eY + step * k;
    euler(:, i) = eY(1:2);

    % Heun
    k1 = f(0, hY);
    k2 = f(0, hY + step * k1);
    hY = hY + (step / 2) * (k1 + k2);
    heun(:, i) = hY(1:2);
end

% Plot trajectory
figure()
plot(heun(1,:), heun(2,:), ...
     euler(1,:), euler(2,:), '--', ...
     oy(:,1), oy(:,2), '-.')
legend({'Heun', 'Euler', 'ode45'}, 'Location', 'best')
xlabel('x')
ylabel('y')
title('Trajectory for pendulum, \Deltat = 0.00001')
snapnow

%% Comment on plots (c)
% We can see from the plots for $\Delta t = 0.0001$ and $\Delta t =
% 0.00001$ that Euler's method seems to improve by a lot - converging toward
% the other curves, while we cannot see any real change in Heun's method, which is
% already good enough.

%% Time evolution of mechanical energy (e)
% Set new time step and reset initial conditions
step = 0.001;
N = ceil(endt / step);
step = endt / N;

eY = [x0; y0; u0; v0];
eE = zeros(1, N);
eE(1) = -((x0 - R)^2 + y0^2)^(-1/2) - ((x0 + R)^2 + y0^2)^(-1/2);

hY = [x0; y0; u0; v0];
hE = eE;

t = linspace(0, 10, N);

% Compute mechanical energy evolution
for i=2:N
    % Euler
    k = f(0, eY);
    eY = eY + step * k;
    eE(i) = - 1 / sqrt((eY(1) - R)^2 + eY(2)^2 + d^2) ...
            - 1 / sqrt((eY(1) + R)^2 + eY(2)^2 + d^2) ...
            + (eY(3)^2 + eY(4)^2) / 2;

    % Heun
    k1 = f(0, hY);
    k2 = f(0, hY + step * k1);
    hY = hY + (step / 2) * (k1 + k2);
    hE(i) = - 1 / sqrt((hY(1) - R)^2 + hY(2)^2 + d^2) ...
            - 1 / sqrt((hY(1) + R)^2 + hY(2)^2 + d^2) ...
            + (hY(3)^2 + hY(4)^2) / 2;
end

% Plot mechanical energy versus time
figure()
plot(t, eE, t, hE, '--')
legend({'Euler', 'Heun'}, 'Location', 'best')
xlabel('t')
ylabel('Mechanical energy')
title('Total mechanical energy in Euler and Heun schemes')
snapnow

%% Comment on plot (e)
% We see from the figure that the Euler scheme does not conserve energy.
% Probably the Heun scheme does not either, but it is much more stable at
% the very least.

%% Chaotic nature (g)

% Set initial values
x0 = 2.0;
y01 = 0.1;
y02 = 0.10001;
u0 = 0.1;
v0 = 0.0;

[t1, y1] = ode45(@f, [0 10], [x0, y01, u0, v0]);
[t2, y2] = ode45(@f, [0 10], [x0, y02, u0, v0]);

figure()
plot(y1(:,1), y1(:,2), y2(:,1), y2(:,2), '--')
legend({'y0 = 0.1', 'y0 = 0.10001'}, 'Location', 'best')
xlabel('x')
ylabel('y')
title('Erratic behaviour upon small change in ICs')
snapnow

%% RHS of ODE (function f)
function Y = f(t, y)
global R d
    Y = zeros(4, 1);
    Y(1) = y(3);
    Y(2) = y(4);
    Y(3) = - (y(1) - R) * ((y(1) - R)^2 + y(2)^2 + d^2)^(-3/2) ...
           - (y(1) + R) * ((y(1) + R)^2 + y(2)^2 + d^2)^(-3/2);
    Y(4) = - y(2) * ((y(1) - R)^2 + y(2)^2 + d^2)^(-3/2) ...
           - y(2) * ((y(1) + R)^2 + y(2)^2 + d^2)^(-3/2);
    
end

