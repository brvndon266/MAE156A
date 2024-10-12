%% Pendulum Estimate Swing Time using Runge-Kutta Method 

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

% Function to calculate angular acceleration (alpha)
alpha_func = @(theta, omega, r) (-m * g * L * sin(theta) - mu * sign(omega) * m * g * r) / (m * L^2);

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

    % Runge-Kutta 4th order method loop
    for n = 1:length(t)-1
        % Adjust friction for specific radii
        if r == 0.125 * 0.0254  % for 1/8 inch radius
            friction_factor = 1.0;  % no additional frication adjustment
        elseif r == 0.25 * 0.0254  % for 1/4 inch radius
            friction_factor = 1.0;  % slight adjustment
        else
            friction_factor = 1.0;  % default
        end
        
        % Apply the friction adjustment factor
        alpha_func_adjusted = @(theta, omega) alpha_func(theta, omega, r) * friction_factor;
        
        % Runge-Kutta 4th order calculations
        k1_theta = omega(n);
        k1_omega = alpha_func_adjusted(theta(n), omega(n));
        
        k2_theta = omega(n) + 0.5 * dt * k1_omega;
        k2_omega = alpha_func_adjusted(theta(n) + 0.5 * dt * k1_theta, omega(n) + 0.5 * dt * k1_omega);
        
        k3_theta = omega(n) + 0.5 * dt * k2_omega;
        k3_omega = alpha_func_adjusted(theta(n) + 0.5 * dt * k2_theta, omega(n) + 0.5 * dt * k2_omega);
        
        k4_theta = omega(n) + dt * k3_omega;
        k4_omega = alpha_func_adjusted(theta(n) + dt * k3_theta, omega(n) + dt * k3_omega);
        
        theta(n+1) = theta(n) + (dt / 6) * (k1_theta + 2 * k2_theta + 2 * k3_theta + k4_theta);
        omega(n+1) = omega(n) + (dt / 6) * (k1_omega + 2 * k2_omega + 2 * k3_omega + k4_omega);
        
        % Calculate potential and kinetic energy
        height = L * (1 - cos(theta(n+1)));  % height above the lowest point
        PE(n+1) = m * g * height;  % Potential Energy
        KE(n+1) = 0.5 * m * (L * omega(n+1))^2;  % Kinetic Energy
        TotalEnergy(n+1) = PE(n+1) + KE(n+1);  % Total Energy
        
        % Stop if the pendulum comes to rest
        if abs(omega(n+1)) < 5e-5 && abs(theta(n+1)) < 0.0025
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
