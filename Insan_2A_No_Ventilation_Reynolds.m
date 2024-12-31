% ---------------------------
% By: Muhammad Insan Kamil Ghifari
% Team members: Pandu Kristian Prayoga Simamora, Abdelrahman Ashraf Gomaa , Amir Faisal bin Shaiful Azuar
% Description:
%   This script processes multiple measurement files (PostProcessed_MeasX.csv),
%   calculates the Reynolds number for each case based on the flow velocity and RPM,
%   and plots various relationships between the Reynolds number, RPM, torque coefficient (K_t),
%   power coefficient (K_q), and efficiency (η). It generates subplots for each relationship 
%   and compares the experimental data across different flow velocities (0, 0.20, 0.35, and 0.70 m/s).
% ---------------------------


clc;
clear;
close all;

% Flow velocity for each case
flowVelocities = [0, 0.20, 0.35, 0.70];  % m/s

% Define constants
c = 0.05; % chord length in meters
r = 0.06; % span of the chord in meters
rho = 1000; 
R = 0.0825; 
D = 2 * R;  
nu = 1.0e-6; % kinematic viscosity in m^2/s (typical value for water)

fileNames = {'PostProcessed_Meas1.csv', 'PostProcessed_Meas2.csv', 'PostProcessed_Meas3.csv', 'PostProcessed_Meas4.csv'};

% Prepare figure for multiple scatter plots
figure;
hold on;

% Use colormap 'jet' for coloring
colormap jet;
colors = colormap;  % Get the colors from the colormap

% Calculate the color index for each case
numCases = length(fileNames);
colorStep = floor(size(colors, 1) / numCases);

% Create subplots for each plot
subplot(2, 2, 1); % Plot 1: Reynolds number vs RPM
hold on;

subplot(2, 2, 2); % Plot 2: Kt vs Reynolds
hold on;

subplot(2, 2, 3); % Plot 3: Kq vs Reynolds
hold on;

subplot(2, 2, 4); % Plot 4: Eta vs Reynolds
hold on;

% Loop through each case and process the data
for caseIdx = 1:numCases
    data = readtable(fileNames{caseIdx}, 'Delimiter', '\t');
    data = data(data.avg_velocity_rpm > -10, :);  % Keep only rows where RPM > -10 (or adjust condition)
    
    % Get the RPM data
    rpm = data.avg_velocity_rpm;  % RPM data from the file
    
    % Calculate omega_n (angular velocity in rad/s) from RPM
    omega_n = rpm * (2 * pi / 60);  % Convert RPM to rad/s
    n = rpm * (1/60); % Convert RPM to RPS

    % Get the flow velocity for the current case
    V = flowVelocities(caseIdx);  % m/s
    
    % Calculate Reynolds number for each RPM using the given formula
    Re_c = (sqrt(V.^2 + (omega_n * r).^2) * c) / nu;
    
    % Calculate efficiency (eta) for each time window in the data
    Kt = -data.avg_force ./ (rho * n.^2 * D.^4); 
    Kq = -data.avg_torque ./ (100 * rho * n.^2 * D.^5); 
    J = V ./ (n * D);  % Advance coefficient
    eta = (Kt .* J) ./ (Kq * 2 * pi);

    % Plot Reynolds number vs RPM (subplot 1)
    subplot(2, 2, 1);
    scatter(rpm, Re_c, 36, 'filled', 'DisplayName', ['(V = ', num2str(V), ' m/s)'], 'MarkerFaceColor', colors((caseIdx-1)*colorStep + 1, :));
    xlabel('RPM');
    ylabel('Reynolds Number ($Re$)', 'Interpreter', 'latex');
    title('Reynolds Number vs RPM', 'Interpreter', 'latex');
    grid on;

    % Plot Kt vs Reynolds (subplot 2)
    subplot(2, 2, 2);
    scatter(Re_c, Kt, 36, 'filled', 'DisplayName', ['(V = ', num2str(V), ' m/s)'], 'MarkerFaceColor', colors((caseIdx-1)*colorStep + 1, :));
    xlabel('Reynolds Number ($Re$)', 'Interpreter', 'latex');
    ylabel('$K_t$', 'Interpreter', 'latex');
    title('$K_t$ vs Reynolds Number', 'Interpreter', 'latex');
    grid on;
    ylim([0, 0.3]); % Set y-axis to start from 0
    
    % Plot Kq vs Reynolds (subplot 3)
    subplot(2, 2, 3);
    scatter(Re_c, Kq, 36, 'filled', 'DisplayName', ['(V = ', num2str(V), ' m/s)'], 'MarkerFaceColor', colors((caseIdx-1)*colorStep + 1, :));
    xlabel('Reynolds Number ($Re$)', 'Interpreter', 'latex');
    ylabel('$K_q$', 'Interpreter', 'latex');
    title('$K_q$ vs Reynolds Number', 'Interpreter', 'latex');
    grid on;
    ylim([0, 0.04]); % Set y-axis to start from 0
    
    % Plot Eta vs Reynolds (subplot 4)
    subplot(2, 2, 4);
    scatter(Re_c, eta, 36, 'filled', 'DisplayName', ['(V = ', num2str(V), ' m/s)'], 'MarkerFaceColor', colors((caseIdx-1)*colorStep + 1, :));
    xlabel('Reynolds Number ($Re$)', 'Interpreter', 'latex');
    ylabel('Efficiency ($\eta$)', 'Interpreter', 'latex');
    title('Efficiency vs Reynolds Number', 'Interpreter', 'latex');
    grid on;
    ylim([0, 1]); % Set y-axis to start from 0
end

% Set legend
legend('show', 'Interpreter', 'latex', 'Location', 'best');

% Hold off to stop further plotting
hold off;