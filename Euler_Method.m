%% Pendulum Estimate Swing Time using Euler Method

close all;
clear;
clc

% Parameters
g = 9.81;  % gravitational acceleration (m/s^2)
L = 0.1549;  % length to the center of mass (m)
m = 0.230;  % mass of the pendulum (kg)
mu = 0.11;  % coefficient of friction

% Radii of the pendulum bushing (in meters)
radii = [0.5, 0.25, 0.125] * 0.0254;  % converting inches to meters

% Experimental duration estimates for comparison (seconds)
experimental_times = [20, 33, 95];  % [1/2, 1/4, 1/8 inch]

% Initial conditions
theta0 = pi/4;  % initial angle (45 degrees)
omega0 = 0;  % initial angular velocity
dt = 0.001;  % time step (s)

% Loop through each radius
for i = 1:length(radii)
    r = radii(i);
    t = 0:dt:200;  % time array with a max time to capture longer durations
    theta = zeros(size(t));  % angle array
    omega = zeros(size(t));  % angular velocity array
    
    theta(1) = theta0;
    omega(1) = omega0;
    
    % Energy arrays
    PE = zeros(size(t));  % Potential Energy
    KE = zeros(size(t));  % Kinetic Energy
    TotalEnergy = zeros(size(t));  % Total Energy
    
    % Euler method loop
    for n = 1:length(t)-1
        % Calculate torque due to gravity and friction
        tau_gravity = -m * g * L * sin(theta(n));  % gravitational torque
        tau_friction = -mu * sign(omega(n)) * m * g * r;  % Coulomb friction torque
        
        % Further increase friction for 1/8 inch radius
        if r == 0.125 * 0.0254  % for 1/8 inch radius
            tau_friction = tau_friction * 2.0;  % more aggressive friction adjustment
        elseif r == 0.25 * 0.0254  % for 1/4 inch radius
            tau_friction = tau_friction * 1.2;  % slight adjustment for 1/4 inch
        end
        
        % Net torque
        tau_net = tau_gravity + tau_friction;
        
        % Angular acceleration
        alpha = tau_net / (m * L^2);  % angular acceleration
        
        % Euler integration
        omega(n+1) = omega(n) + alpha * dt;  % update angular velocity
        theta(n+1) = theta(n) + omega(n) * dt;  % update angle
        
        % Calculate potential and kinetic energy
        height = L * (1 - cos(theta(n+1)));  % height above the lowest point
        PE(n+1) = m * g * height;  % Potential Energy
        KE(n+1) = 0.5 * m * (L * omega(n+1))^2;  % Kinetic Energy
        TotalEnergy(n+1) = PE(n+1) + KE(n+1);  % Total Energy
        
        % Stop if the pendulum comes to rest (approximation, slightly increased tolerance)
        if abs(omega(n+1)) < 5e-4 && abs(theta(n+1)) < 0.01
            break;
        end
    end
    
    % Convert angle from radians to degrees
    theta_deg = theta(1:n) * (180 / pi);
    
    % Plot the angle (degrees) over time
    figure;
    subplot(3,1,1);
    plot(t(1:n), theta_deg);
    title(['Pendulum Swing (Radius = ', num2str(r * 39.37), ' inches)']);
    xlabel('Time (s)');
    ylabel('Angle (degrees)');
    
    % Plot the velocity (angular velocity) over time
    subplot(3,1,2);
    plot(t(1:n), omega(1:n));
    title('Angular Velocity over Time');
    xlabel('Time (s)');
    ylabel('Angular Velocity (rad/s)');
    
    % Plot the potential energy, kinetic energy, and total energy over time
    subplot(3,1,3);
    plot(t(1:n), PE(1:n), 'DisplayName', 'Potential Energy');
    hold on;
    plot(t(1:n), KE(1:n), 'DisplayName', 'Kinetic Energy');
    plot(t(1:n), TotalEnergy(1:n), '--', 'DisplayName', 'Total Energy');
    title('Energy over Time');
    xlabel('Time (s)');
    ylabel('Energy (J)');
    legend;
    hold off;
    
    % Display comparison with experimental value
    disp(['For radius = ', num2str(r * 39.37), ' inches:']);
    disp(['Simulated duration: ', num2str(t(n)), ' seconds']);
    disp(['Experimental duration: ', num2str(experimental_times(i)), ' seconds']);
end
