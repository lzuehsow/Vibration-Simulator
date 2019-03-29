function suspension()
m = 4000; % kg
r2 = 0.64; % m
k1 = 20000; % N/m
k2 = 20000; % N/m
l1 = 0.9; % m
l2 = 1.4; % m

I = m * r2;
M = [m 0; 0 I];
K = [k1+k2 (-k1*l1 + k2*l2); (-k1*l1 + k2*l2) ((l1^2)*k1 + (l2^2)*k2)];
C = [4000 1000; 1000 5540];

L = chol(M); % m ^ 0.5
L_inv = inv(L);

K_tilda = L_inv * K * L_inv;
C_tilda = L_inv * C * L_inv;
[eigvec,eigval] = eig(K_tilda);

P = eigvec;

diag = transpose(P) * K_tilda * P;
[V,D] = eig(diag);

%% ODE

% Define simulation parameters
t_span = linspace(0,10);
x_init = [0.2 0; 0.1 0];

x0 = x_init(1); % m
theta0 = x_init(2); % rad
vx0 = x_init(3); % m/s
vtheta0 = x_init(4); % rad/s

Z_0 = [x0, theta0, vx0, vtheta0];     %specify initial conditions (rads, rads/sec)

[Time, zout] = ode45(@ode, t_span, Z_0);

function res = ode(~, Z)
res = zeros(4,1);
res(1) = Z(3); % dx2dt
res(2) = Z(4); % dTheta2dt

% % Undamped
% res(3) = (-(k1 + k2) * Z(1) + (k1*l1 - k2*l2) * Z(2)) / m;
% res(4) = ((l1*k1 - l2*k2) * Z(1) - ((l1^2)*k1 + (l2^2)*k2) * Z(2)) / I;

K = [k1+k2 (-k1*l1 + k2*l2); (-k1*l1 + k2*l2) ((l1^2)*k1 + (l2^2)*k2)];
res(3) = (-K(1) * Z(1) + -K(2) * Z(2)) / m;
res(4) = (-K(3) * Z(1) - K(4) * Z(2)) / I;
end

%% Analytical Soln
% Initial conditions in x and q space
q_init = transpose(P) * L * x_init;

% Find damping constant
zeta = ((C / K) * eigval) / (2 * sqrt(eigval));

% Find damped natural frequency
Wd = sqrt(eigval) * sqrt ([1 0; 0 1] - (zeta.^2));

% Equation of motion in q space
undamped = [q_init(1) * cos(sqrt(eigval(1)) * t_span) + (q_init(3)/sqrt(eigval(1))) * sin(sqrt(eigval(2)) * t_span); ...
q_init(2) * cos(sqrt(eigval(4)) * t_span) + (q_init(4)/sqrt(eigval(4))) * sin(sqrt(eigval(2)) * t_span)];

underdamped = [exp(-zeta(1) * sqrt(eigval(1)) * t_span) .* ...
    (q_init(1) * cos(Wd(1) * t_span) + ((q_init(2) + zeta(1) * sqrt(eigval(1)) * q_init(1)) / Wd(1)) * sin(Wd(1) * t_span)); ...
    exp(-zeta(4) * sqrt(eigval(4)) * t_span) .* ...
    (q_init(3) * cos(Wd(4) * t_span) + ((q_init(4) + zeta(4) * sqrt(eigval(3)) * q_init(3)) / Wd(4)) * sin(Wd(4) * t_span))];

% Equation of motion converted to x space
undamped_analytical = L_inv * P * undamped;
damped_analytical = L_inv * P * underdamped;

%% Plotting

figure;
hold on;

plot(Time, zout(:, 1, 1, 1), Time, damped_analytical(1,:));

% plot(Time, zout(:, 1, 1, 1), Time, zout(:, 2, 1, 1)); % ODE
% legend('Theta', 'X')

plot(t_span, damped_analytical);               % Analytical Soln
legend('X', 'Theta', 'X Analytical', 'Theta Analytical');

title('X & Theta vs Time');
xlabel('Time (s)');
end