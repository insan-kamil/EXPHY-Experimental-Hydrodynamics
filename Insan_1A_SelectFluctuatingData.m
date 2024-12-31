% ---------------------------
% By: Muhammad Insan Kamil Ghifari
% Team members: Pandu Kristian Prayoga Simamora, Abdelrahman Ashraf Gomaa , Amir Faisal bin Shaiful Azuar
% Description: 
%   This script allows the user to interactively select a raw measurement file 
%   (MeasX.csv), plot key data columns from the selected file, and then select 
%   time windows of fluctuating data using mouse clicks. The script extracts raw 
%   data from the selected time windows and calculates the mean RPM and average 
%   water speed for each window. The results are saved to both the base workspace 
%   and a CSV file (FluctuatingData_MeasX.csv).
% ---------------------------

clc;
clear;
close all;

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

% Step 5: Extract raw data for each selected time window
for i = 1:size(windows, 1)
    time_window_start = windows(i, 1);
    time_window_end = windows(i, 2);
    
    % Extract the raw data within the selected window
    steady_data = data(data.Time_ms >= time_window_start & data.Time_ms <= time_window_end, :);
    
    % Add a new column for the time window index
    window_index = repmat(i, height(steady_data), 1);  % Create a column with the current window index
    
    % Add the window index as a new column
    steady_data.Window_Index = window_index;
    
    % Keep only the necessary columns: Time_ms, Velocity_Rpm, WaterSpeed1_m_s, Force_N, Torque_N_cm, Window_Index
    steady_data = steady_data(:, {'Time_ms', 'Velocity_Rpm', 'WaterSpeed1_m_s', 'Force_N', 'Torque_N_cm', 'Window_Index'});
    
    % Append the current window data to the result table
    result_table = [result_table; steady_data];
end

% Step 6: Create the new column for mean RPM
result_table.Mean_Velocity_Rpm = NaN(height(result_table), 1);  % Initialize the column with NaN values

% Step 7: Calculate the mean RPM for each window index and store it in the new column
% Group by Window_Index and calculate the mean RPM for each window
mean_rpm_per_window = varfun(@mean, result_table, 'InputVariables', 'Velocity_Rpm', 'GroupingVariables', 'Window_Index');

% Step 8: Update the `Mean_Velocity_Rpm` column with the calculated mean values
for i = 1:height(mean_rpm_per_window)
    % Get the mean RPM and round it to the nearest 5
    rounded_mean_rpm = round(mean_rpm_per_window.mean_Velocity_Rpm(i) / 5) * 5;
    
    % Assign the rounded value to the corresponding rows in the result table
    result_table(result_table.Window_Index == mean_rpm_per_window.Window_Index(i), 'Mean_Velocity_Rpm') = ...
        {rounded_mean_rpm};  % Use {} to assign scalar value
end

% Step 9: Calculate average water speed (WaterSpeed1_m_s) for each time window
% Group by Window_Index and calculate the mean of WaterSpeed1_m_s for each window
mean_water_speed_per_window = varfun(@mean, result_table, 'InputVariables', 'WaterSpeed1_m_s', 'GroupingVariables', 'Window_Index');

% Step 10: Add a new column for Average Water Speed (WaterSpeed1_m_s)
result_table.Avg_WaterSpeed1_m_s = NaN(height(result_table), 1);  % Initialize with NaN values

% Step 11: Update the `Avg_WaterSpeed1_m_s` column with the calculated averages
for i = 1:height(mean_water_speed_per_window)
    % Get the average water speed and assign it to the table (using {})
    result_table(result_table.Window_Index == mean_water_speed_per_window.Window_Index(i), 'Avg_WaterSpeed1_m_s') = ...
        {mean_water_speed_per_window.mean_WaterSpeed1_m_s(i)};  % Use {} to assign scalar value
end

% Step 12: Save the updated table to a CSV file
csv_filename = sprintf('FluctuatingData_Meas%d.csv', selected_case);  % FluctuatingData_MeasX.csv
writetable(result_table, csv_filename, 'Delimiter', '\t');
