% ---------------------------
% By: Muhammad Insan Kamil Ghifari
% Team members: Pandu Kristian Prayoga Simamora, Abdelrahman Ashraf Gomaa , Amir Faisal bin Shaiful Azuar
% Description:
%   This script processes multiple measurement files (FluctuatingData_MeasX.csv) 
%   related to propeller performance under varying ventilation conditions, specifically 
%   focused on different values of `h/r` (distance between propeller and hull radius). 
%   The script performs the following tasks:
%   1. It calculates key performance parameters such as Thrust, Torque, and Efficiency (η) 
%      by processing the raw data, subtracting offsets, and calculating the advance coefficient (J).
%   2. The script generates boxplots for Thrust, Torque, and Efficiency as a function of `J`, 
%      comparing the results for different ventilation cases.
%   3. A table is generated for each case, which includes statistics such as mean, median, 
%      and standard deviation (STD) of Thrust and Torque at each RPM.
%   4. The tables are saved to CSV files for further analysis.
% ---------------------------

%% Plot
clc;
clear;
close all;

% Case names, colors, and file names
caseNames = {'h/r = 2.03', 'h/r = 1.24', 'h/r = 1.36'};
colors = [0, 1, 0; 1, 0, 0; 0, 0, 1];  % Green, Red, Blue
file_names = {'FluctuatingData_Meas3.csv', 'FluctuatingData_Meas5.csv', 'FluctuatingData_Meas6.csv'};

% Constants and offsets
R = 0.0825;  % Radius in meters
D = 2 * R;   % Diameter in meters
reference_water_speed = 0.35;  % Average water speed (m/s)
force_offset = 1.80879894545455;
torque_offset = 1.63032294641148;

% Create figures for Thrust and Torque boxplots
for fig_idx = 1:2
    figure;
    for i = 1:length(file_names)
        % Read and process data
        data = readtable(file_names{i}, 'Delimiter', '\t', 'VariableNamingRule', 'modify');
        data = data(data.Mean_Velocity_Rpm > 110, :);  % Keep only rows with Mean_Velocity_Rpm > 110
        data.Force_N = data.Force_N - force_offset;  % Subtract force offset
        data.Torque_N_cm = data.Torque_N_cm - torque_offset;  % Subtract torque offset
        
        % Calculate 'J'
        n = data.Mean_Velocity_Rpm * (1 / 60);  % RPM to Hz
        J = reference_water_speed ./ (n * D);  % Calculate J
        
        % Determine whether to plot Thrust or Torque
        if fig_idx == 1
            % Thrust: Plot Thrust (T) as a function of J
            values = -data.Force_N;
            ylabel_text = '$T$ (N)';
            title_text = ['$\textbf{Thrust vs J}$ - ' caseNames{i}];
            % Set y-axis limits for Thrust
            ylim_values = [0, 40];
        else
            % Torque: Plot Torque (Q) as a function of J
            values = -data.Torque_N_cm / 100;  % Convert from N.cm to N.m
            ylabel_text = '$Q$ (N.m)';
            title_text = ['$\textbf{Torque vs J}$ - ' caseNames{i}];
            % Set y-axis limits for Torque
            ylim_values = [0, 0.75];
        end
        
        % Plot the boxplot for the current case
        subplot(1, 3, i);
        h = boxplot(values, J, 'Whisker', Inf);
        set(h, 'Color', colors(i, :));  % Set color for all elements (whiskers, box, median)
        
        % Format the plot
        title(title_text, 'Interpreter', 'latex');
        xlabel('$J$ (dimensionless)', 'Interpreter', 'latex');
        ylabel(ylabel_text, 'Interpreter', 'latex');
        grid on;
        
        % Set y-axis limits
        ylim(ylim_values);
        
    end
end

%% Table

% Create a cell array to store tables for each case
tables = cell(1, 3);  % Create a cell to store tables for each case

% Loop through each case to process the data
for i = 1:length(file_names)
    % Read and process data
    data = readtable(file_names{i}, 'Delimiter', '\t', 'VariableNamingRule', 'modify');
    data = data(data.Mean_Velocity_Rpm > 110, :);  % Keep only rows with Mean_Velocity_Rpm > 110
    
    data.Force_N = -1*(data.Force_N - force_offset);  % Subtract force offset
    data.Torque_N_cm = -1*(data.Torque_N_cm - torque_offset);  % Subtract torque offset
    
    % Calculate 'J'
    n = data.Mean_Velocity_Rpm * (1 / 60);  % RPM to Hz
    J = reference_water_speed ./ (n * D);  % Calculate J
    
    % Initialize variables for statistics
    mean_rpm = unique(data.Mean_Velocity_Rpm);  % Get unique RPM values
    num_rpm = length(mean_rpm);  % Number of unique RPM values
    
    % Initialize the table for this case
    case_table = table('Size', [num_rpm, 7], 'VariableTypes', {'double', 'double', 'double', 'double', 'double', 'double', 'double'}, ...
        'VariableNames', {'Mean_RPM', 'Mean_Thrust', 'Median_Thrust', 'STD_Thrust', 'Mean_Torque', 'Median_Torque', 'STD_Torque'});
    
    % Loop through each unique RPM
    for j = 1:num_rpm
        % Extract data for the current RPM
        rpm_data = data(data.Mean_Velocity_Rpm == mean_rpm(j), :);
        
        % Calculate mean, median, and standard deviation for Thrust and Torque
        mean_thrust = mean(rpm_data.Force_N);
        median_thrust = median(rpm_data.Force_N);
        std_thrust = std(rpm_data.Force_N);
        
        mean_torque = mean(rpm_data.Torque_N_cm / 100);  % Convert Torque from N.cm to N.m
        median_torque = median(rpm_data.Torque_N_cm / 100);  % Convert Torque from N.cm to N.m
        std_torque = std(rpm_data.Torque_N_cm / 100);  % Convert Torque from N.cm to N.m
        
        % Populate the table for this RPM
        case_table(j, :) = {mean_rpm(j), mean_thrust, median_thrust, std_thrust, mean_torque, median_torque, std_torque};
    end
    
    % Store the table for the current case
    tables{i} = case_table;
end

% Define the directory where CSV files will be saved
output_dir = 'OutputTables';  % You can change this to any valid path

% Create the directory if it doesn't exist
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% Save the tables to CSV files
for i = 1:length(caseNames)
    % Clean up the case name to avoid issues with file naming
    cleaned_case_name = strrep(caseNames{i}, '/', '_');  % Replace '/' with '_'
    cleaned_case_name = strrep(cleaned_case_name, ' ', '_');  % Replace spaces with '_'
    
    % Generate a filename for each case (save inside the output directory)
    csv_filename = fullfile(output_dir, sprintf('VariationData_%s.csv', cleaned_case_name));
    
    % Save the table to a CSV file
    writetable(tables{i}, csv_filename);
    
    % Optionally, display a message confirming the save
    disp(['Table for ', caseNames{i}, ' saved as ', csv_filename]);
end
