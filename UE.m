num_known_interferers = 5;
num_unknown_interferers = 40;

% Known Interferers - Fixed Locations
radii = [60, 100, 140, 180, 220]; % Fixed distances
angles = 2*pi*rand(1, num_known_interferers); % Random angles
x_known = radii .* cos(angles);
y_known = radii .* sin(angles);

% Unknown Interferers - Random Locations
r_unknown = 250 + (500-250) * rand(1, num_unknown_interferers);
theta_unknown = 2 * pi * rand(1, num_unknown_interferers);
x_unknown = r_unknown .* cos(theta_unknown);
y_unknown = r_unknown .* sin(theta_unknown);

% Plot
figure; hold on; grid on;
scatter(0, 0, 100, 'ro', 'filled'); % Serving BS
scatter(100, 0, 100, 'bo', 'filled'); % Desired UE
scatter(x_known, y_known, 100, 'gs', 'filled'); % Known Interferers
scatter(x_unknown, y_unknown, 50, 'k*'); % Unknown Interferers
legend('Serving BS', 'Desired UE', 'Known Interferers', 'Unknown Interferers');
title('UE and Interferer Placement');
xlabel('X (meters)');
ylabel('Y (meters)');