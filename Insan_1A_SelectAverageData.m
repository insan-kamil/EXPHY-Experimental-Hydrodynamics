% ---------------------------
% By: Muhammad Insan Kamil Ghifari
% Team members: Pandu Kristian Prayoga Simamora, Abdelrahman Ashraf Gomaa , Amir Faisal bin Shaiful Azuar
% Description: 
%   This script allows the user to interactively select a raw measurement file 
%   (MeasX.csv), plot the key data columns from the selected file, and then 
%   select time windows of steady-state data using mouse clicks. The script calculates 
%   the average values of Water Speed, Motor Velocity, Force, and Torque for each 
%   selected time window and saves the results to both the base workspace and a CSV file 
%   (PostProcessed_MeasX.csv) for further analysis.
% ---------------------------

close all;
clc;
clear;

file_names = {'Meas1.csv', 'Meas2.csv', 'Meas3.csv', 'Meas4.csv', 'Meas5.csv', 'Meas6.csv'};
column_names = {'Time_ms', 'Velocity_Rpm', 'WaterSpeed1_Volts', 'WaterSpeed2_Volts', 'Force_mV_V', ...
    'Torque_mV_V', 'WaterSpeed1_m_s', 'WaterSpeed2_m_s', 'Force_N', 'Torque_N_cm'};

% Ask the user to choose which case to plot
fprintf('Available cases:\n');
for i = 1:length(file_names)
    fprintf('%d: %s\n', i, file_names{i});
end

selected_case = input('Enter the case number to plot: ');

% Validate the input
if selected_case < 1 || selected_case > length(file_names)
    error('Invalid case number selected. Please choose a valid case number.');
end

% Step 1: Read the selected file
filename = file_names{selected_case};
opts = detectImportOptions(filename);
opts.Delimiter = '\t';
opts.VariableNamingRule = 'modify';
data = readtable(filename, opts);
data.Properties.VariableNames = column_names;  % Assign column names

% Step 2: Plot the selected case (full screen figure)
figure('units', 'normalized', 'outerposition', [0 0 1 1]);  % Full screen figure
set(gcf, 'Name', filename, 'NumberTitle', 'off');  % Set figure window title to the filename

% Rearranged subplots to place Water Speed at the top
subplot(4, 1, 1);
plot(data.Time_ms, data.WaterSpeed1_m_s);
xlabel('Time (ms)');
ylabel('Water Speed 1 (m/s)');
title('Water Speed 1 vs Time');
grid on;

subplot(4, 1, 2);
plot(data.Time_ms, data.Velocity_Rpm);
xlabel('Time (ms)');
ylabel('Motor Velocity (rpm)');
title('Motor Velocity vs Time');
grid on;

subplot(4, 1, 3);
plot(data.Time_ms, data.Force_N);
xlabel('Time (ms)');
ylabel('Force (N)');
title('Force vs Time');
grid on;

subplot(4, 1, 4);  % Added subplot for Torque
plot(data.Time_ms, data.Torque_N_cm);
xlabel('Time (ms)');
ylabel('Torque (N.cm)');
title('Torque vs Time');
grid on;

% Step 3: Allow the user to select multiple steady-state windows using ginput
disp('Select steady-state windows by clicking two points for each window (start and end).');
disp('Press Enter when done.');

windows = [];  % To store the start and end points of each window
while true
    [x1, ~] = ginput(1);  % Let user select the first point for the window
    if isempty(x1)  % Break if the user presses Enter (no more clicks)
        break;
    end
    
    % Draw vertical line for the first click (show immediately after the first click, in green)
    subplot(4, 1, 1);
    line([x1 x1], ylim, 'Color', 'g', 'LineStyle', '--');
    subplot(4, 1, 2);
    line([x1 x1], ylim, 'Color', 'g', 'LineStyle', '--');
    subplot(4, 1, 3);
    line([x1 x1], ylim, 'Color', 'g', 'LineStyle', '--');
    subplot(4, 1, 4);
    line([x1 x1], ylim, 'Color', 'g', 'LineStyle', '--');
    
    % Wait for the second click (end of the time window)
    [x2, ~] = ginput(1);  % Let user select the second point for the window
    
    % Draw vertical line for the second click (immediately after second click, in red)
    subplot(4, 1, 1);
    line([x2 x2], ylim, 'Color', 'r', 'LineStyle', '--');
    subplot(4, 1, 2);
    line([x2 x2], ylim, 'Color', 'r', 'LineStyle', '--');
    subplot(4, 1, 3);
    line([x2 x2], ylim, 'Color', 'r', 'LineStyle', '--');
    subplot(4, 1, 4);
    line([x2 x2], ylim, 'Color', 'r', 'LineStyle', '--');
    
    % Store the selected window (start and end points)
    windows = [windows; min(x1, x2), max(x1, x2)];
end

% Step 4: Create a table to store the results
result_table = [];  % Table to store the results

% Step 5: Calculate and store the averages for each selected time window
for i = 1:size(windows, 1)
    time_window_start = windows(i, 1);
    time_window_end = windows(i, 2);
    
    % Extract the data within the selected window
    steady_data = data(data.Time_ms >= time_window_start & data.Time_ms <= time_window_end, :);
    
    % Calculate the averages for the selected window
    avg_water_speed = mean(steady_data.WaterSpeed1_m_s);
    avg_velocity_rpm = mean(steady_data.Velocity_Rpm);
    avg_force = mean(steady_data.Force_N);
    avg_torque = mean(steady_data.Torque_N_cm);
    
    % Append the results for the current window to the table
    result_table = [result_table; table(time_window_start, time_window_end, avg_water_speed, avg_velocity_rpm, avg_force, avg_torque)];
end

% Step 6: Name the table according to the selected Meas case
table_name = sprintf('%s_results', file_names{selected_case}(1:end-4));  % Remove '.csv' and add '_results'
assignin('base', table_name, result_table);  % Assign the table to the base workspace with the case name

% Step 7: Save the table to a CSV file named 'PostProcessed_<MeasX>.csv'
csv_filename = sprintf('PostProcessed_%s.csv', file_names{selected_case}(1:end-4));  % PostProcessed_Meas1.csv

% Write the table to CSV with tab separation
writetable(result_table, csv_filename, 'Delimiter', '\t');

% Display the table in the command window
disp('Results for the selected time windows:');
disp(result_table);
