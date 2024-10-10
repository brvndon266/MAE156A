% Pendulum Estimate Swing Time using Euler Method

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
        
        % Stop if the pendulum comes to rest (approximation, slightly increased tolerance)
        if abs(omega(n+1)) < 5e-4 && abs(theta(n+1)) < 0.01
            break;
        end
    end
    
    % Calculate total duration of swing
    total_time = t(n);
    
    % Plot the results for each radius
    figure;
    plot(t(1:n), theta(1:n));
    title(['Pendulum Swing with Radius = ', num2str(r * 39.37), ' inches']);
    xlabel('Time (s)');
    ylabel('Angle (rad)');
    
    % Display comparison with experimental value
    disp(['For radius = ', num2str(r * 39.37), ' inches:']);
    disp(['Simulated duration: ', num2str(total_time), ' seconds']);
    disp(['Experimental duration: ', num2str(experimental_times(i)), ' seconds']);
end
