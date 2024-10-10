%% Pendulum Estimate Swing Time using Euler Method

close all;
clear;
clc;

% Initialize Parameters
theta_init=deg2rad(45); % units: [radians]
thetadot_init=0;        % units: [radians/second]

g=9.81; % gravity (m/s^2)
L=.167;    % Length from center of mass to pivot point (m)
theta=theta_init;
thetadot=thetadot_init;


% Store the intial values of the angle, angular velocity and the time
theta_Euler(1)=theta;
thetadot_Euler(1)=thetadot;
t_Euler(1)=0;

tic % Start timer in order to track how long will the pendulum swing

for i=2:point % starts from point 2 until the pendulum stops oscillating
    t= t + t_step;
    newtheta=theta+ thetadot*t_step; % New angle after each iteration
    theta2dot=-(g*L)*sin(theta);     % Angular acceleration
    newthetadot= thetadot + theta2dot*t_step; % New angular velocity
% Update values to new values

    theta=newtheta;
    thetadot=newthetadot;

    t_Euler(i)=t;
    theta_Euler(1)=theta;
    thetadot_Euler(1)=thetadot;
end

toc % Stops time measurement and displays time









