% ---------------------------
% By: Muhammad Insan Kamil Ghifari
% Team members: Pandu Kristian Prayoga Simamora, Abdelrahman Ashraf Gomaa , Amir Faisal bin Shaiful Azuar
% Description:
%   This script processes multiple post-processed measurement files for lower immersion cases (PostProcessed_MeasX.csv) 
%   and calculates key parameters such as torque coefficient (K_t), power coefficient (K_q), and efficiency (η) 
%   for each case. It then compares the experimental data to theoretical reference data, 
%   plotting the results on a series of scatter plots for different cases with distinct colors and marker shapes.
%   The script generates subplots to show the relationships between these parameters and the advance coefficient (J).
% ---------------------------

clc;
clear;
close all;

rho = 1000; 
R = 0.0825; 
D = 2 * R;  % Define the diameter of the propeller

% Define the case names for each file
caseNames = {'h/R = 2.03', 'h/R = 1.24', 'h/R = 1.36'};

% Define the file names
fileNames = {'PostProcessed_Meas3.csv', 'PostProcessed_Meas5.csv', 'PostProcessed_Meas6.csv'};

% Water speed for all cases
avg_water_speed_current_file = 0.35;

% Theoretical Data (Reference)
J_theo = linspace(0, 0.77, 1000);
K_T_theo = 0.2469 - (0.2458 .* J_theo) - (0.2258 .* (J_theo.^2)) + (0.0675 .* (J_theo.^3));
K_Q_theo = 0.02394 - 0.01899 * J_theo - 0.01384 * J_theo.^2 - 0.001973 * J_theo.^3;
eta_theo = (K_T_theo .* J_theo) ./ (K_Q_theo * 2 * pi);

% Colors for cases (using predefined colors)
colors = [
    0, 1, 0;   % Green
    1, 0, 0;   % Red
    0, 0, 1;   % Blue
];

% Loop over each file and process data
% Loop over each file and process data
for i = 1:length(fileNames)
    % Read data
    data = readtable(fileNames{i}, 'Delimiter', '\t');
    data = data(:, 3:end); 
    
    % Filter out rows where avg_velocity_rpm is between 0 and 210
    data = data(data.avg_velocity_rpm < 0 | data.avg_velocity_rpm > 210, :);
    
    % Display the water speed for the current case
    disp(['Water Speed for case ', num2str(i), ' (', caseNames{i}, '): ', num2str(avg_water_speed_current_file)]);
    
    % Calculate parameters
    n = data.avg_velocity_rpm * (1/60); 
    J = avg_water_speed_current_file ./ (n * D); 
    Kt = -data.avg_force ./ (rho * n.^2 * D.^4);
    Kq = -data.avg_torque ./ (100 * rho * n.^2 * D.^5); % cm to m
    eta_exp = (Kt .* J) ./ (Kq * 2 * pi);
    
    markerSize = 50;

    % Assign different marker shapes for each case
    markerShapes = {'o', 's', 'd'}; % Circle, Square, Diamond
    
    % Plot Kt, Kq, and Eta (Efficiency) for each case using the reusable plot function
    plotData(J, Kt, colors(i, :), caseNames{i}, '$K_t$', 1, [0, 0.77], [0, 0.26], markerSize, markerShapes{i});
    plotData(J, Kq, colors(i, :), caseNames{i}, '$K_q$', 2, [0, 0.77], [0, 0.026], markerSize, markerShapes{i});
    plotData(J, eta_exp, colors(i, :), caseNames{i}, '$\eta$', 3, [0, 0.77], [0, 0.6], markerSize, markerShapes{i});
end


% Plot theoretical data after looping through files
plotTheoreticalData(J_theo, K_T_theo, K_Q_theo, eta_theo);

% Function to handle plotting
function plotData(J, data, color, caseName, ylabelText, figNum, xlimVals, ylimVals, markerSize, markerShape)
    figure(figNum);
    hold on;
    scatter(J, data, markerSize, 'MarkerFaceColor', 'none', 'MarkerEdgeColor', color, 'Marker', markerShape, 'DisplayName', ['$' caseName '$']);
    xlabel('$J$ (Advance Coefficient)', 'Interpreter', 'latex');
    ylabel(ylabelText, 'Interpreter', 'latex');
    
    % Correct title formatting using sprintf to ensure correct LaTeX interpretation
    titleStr = sprintf('%s vs $J$', ylabelText);
    title(titleStr, 'Interpreter', 'latex');
    
    xlim(xlimVals);
    ylim(ylimVals);
    legend('show', 'Interpreter', 'latex', 'Location', 'best');
    grid on;
end


% Function to plot theoretical data
function plotTheoreticalData(J_theo, K_T_theo, K_Q_theo, eta_theo)
    % Plot Kt vs J
    figure(1);
    plot(J_theo, K_T_theo, '-', 'LineWidth', 1, 'DisplayName', '$K_t$ (Reference)', 'Color', [0, 0, 0, 1]); % Solid black line
    
    % Plot Kq vs J
    figure(2);
    plot(J_theo, K_Q_theo, '-', 'LineWidth', 1, 'DisplayName', '$K_q$ (Reference)', 'Color', [0, 0, 0, 1]); % Solid black line
    
    % Plot Efficiency vs J
    figure(3);
    plot(J_theo, eta_theo, '-', 'LineWidth', 1, 'DisplayName', '$\eta$ (Reference)', 'Color', [0, 0, 0, 1]); % Solid black line
    
    % Customize all figures
    for figNum = 1:3
        figure(figNum);
        xlabel('$J$ (Advance Coefficient)', 'Interpreter', 'latex');
        legend('show', 'Interpreter', 'latex', 'Location', 'best');
        grid on;
    end
end
