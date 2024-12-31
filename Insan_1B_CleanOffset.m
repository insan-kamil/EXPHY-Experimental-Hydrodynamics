% ---------------------------
% By: Muhammad Insan Kamil Ghifari
% Team members: Pandu Kristian Prayoga Simamora, Abdelrahman Ashraf Gomaa , Amir Faisal bin Shaiful Azuar
% Description: 
%   This script processes multiple post-processed measurement files (PostProcessed_MeasX.csv),
%   applies offsets to the force and torque values based on the first entry of the first file,
%   and calculates and displays the average water speed for each case. 
%   The script updates the force and torque values, saves the updated data back to the file,
%   and generates plots of the average force and torque against the average velocity RPM for each case.
% ---------------------------

clc;
clear;
close all;

fileNames = {'PostProcessed_Meas1.csv', 'PostProcessed_Meas2.csv', 'PostProcessed_Meas3.csv', 'PostProcessed_Meas4.csv', 'PostProcessed_Meas5.csv', 'PostProcessed_Meas6.csv'};

for caseIdx = 1:length(fileNames)
    data = readtable(fileNames{caseIdx}, 'Delimiter', '\t');
    
    % Use the first entry of the first file for the offset values
    if caseIdx == 1
        offsetForce = data.avg_force(1);
        offsetTorque = data.avg_torque(1);
    end
    
    % Apply the offsets to force and torque values
    data.avg_force = data.avg_force - offsetForce;
    data.avg_torque = data.avg_torque - offsetTorque;
    
    avgWaterSpeed = mean(data.avg_water_speed);
    disp(['Case ', num2str(caseIdx), ' (', fileNames{caseIdx}, '): Average Water Speed = ', num2str(avgWaterSpeed)]);
    
    % Write the updated table back to the file
    writetable(data, fileNames{caseIdx}, 'Delimiter', '\t');
    
    % Plot the force and torque data
    figure;
    hold on;
    plot(data.avg_velocity_rpm, data.avg_force, 'o', 'DisplayName', 'Avg Force');
    plot(data.avg_velocity_rpm, data.avg_torque, 'x', 'DisplayName', 'Avg Torque');
    hold off;
    
    title(['Case Identity: ', num2str(avgWaterSpeed), ' m/s']);
    xlabel('Average Velocity RPM');
    ylabel('Force / Torque');
    legend('show');
    grid on;
end
