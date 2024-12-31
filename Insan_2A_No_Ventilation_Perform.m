% ---------------------------
% By: Muhammad Insan Kamil Ghifari
% Team members: Pandu Kristian Prayoga Simamora, Abdelrahman Ashraf Gomaa , Amir Faisal bin Shaiful Azuar
% Description:
%   This script processes four post-processed measurement files for the h/R = 2.03 cases (PostProcessed_MeasX.csv),
%   calculates several performance parameters including the advance coefficient (J), 
%   torque coefficient (K_t), power coefficient (K_q), and efficiency (η) based on the experimental data. 
%   It compares the experimental data with theoretical reference values and generates plots for visual analysis. 
%   The script also fits polynomial curves to the experimental data and calculates the goodness of fit (R²).
% ---------------------------

clc;
clear;
close all;

rho = 1000; 
R = 0.0825; 
D = 2 * R;  
combinedData = table();

for i = 1:4
    fileName = sprintf('PostProcessed_Meas%d.csv', i);
    data = readtable(fileName, 'Delimiter', '\t');
    data = data(:, 3:end); 
    
    avg_water_speed_current_file = mean(data.avg_water_speed);
    % Assign the avg_water_speed_current_file value based on i
    % switch i
    %     case 1
    %         avg_water_speed_current_file = 0;
    %     case 2
    %         avg_water_speed_current_file = 0.20;
    %     case 3
    %         avg_water_speed_current_file = 0.35;
    %     case 4
    %         avg_water_speed_current_file = 0.70;
    % end
    disp(['Water Speed for case ', num2str(i), ': ', num2str(avg_water_speed_current_file)]);
    
    n = data.avg_velocity_rpm * (1/60); 
    J = avg_water_speed_current_file ./ (n * D); 
    Kt = -data.avg_force ./ (rho * n.^2 * D.^4);
    Kq = -data.avg_torque ./ (100 * rho * n.^2 * D.^5); %cm to m
    eta_exp = (Kt .* J) ./ (Kq * 2 * pi);
    
    data.J = J;
    data.Kt = Kt;
    data.Kq = Kq;
    data.eta_exp = eta_exp;
    
    combinedData = [combinedData; data]; 
end

% Remove rows where J is NaN or Inf
cleanData = combinedData(~isnan(combinedData.J) & ~isinf(combinedData.J), :);
% Sort the data based on the J column
sortedData = sortrows(cleanData, 'J');
% Delete the some data
% sortedData(1:3, :) = [];  % First 3 rows removed: Low RPM, still water condition.
sortedData = sortedData(sortedData.avg_velocity_rpm < 0 | sortedData.avg_velocity_rpm > 210, :); % Remove rows where avg_velocity_rpm is between 0 and 310

% Theoretical Data (Reference)
J_theo = linspace(0, 1, 500);
K_T_theo = 0.2469 - (0.2458 .* J_theo) - (0.2258 .* (J_theo.^2)) + (0.0675 .* (J_theo.^3));
K_Q_theo = 0.02394 - 0.01899 * J_theo - 0.01384 * J_theo.^2 - 0.001973 * J_theo.^3;
eta_theo = (K_T_theo .* J_theo) ./ (K_Q_theo * 2 * pi);

% Experimental Data Plot
figure;
ax = gca; % Get the current axes

% Left axis
yyaxis left;
% Reference Kt
plot(J_theo, K_T_theo, '-', 'LineWidth', 2, 'DisplayName', '$K_t$ (Reference)', 'Color', [0, 0.4470, 0.7410, 0.5]); 
hold on;
% Reference Eta
plot(J_theo, eta_theo, '-', 'LineWidth', 2, 'DisplayName', '$\eta$ (Reference)', 'Color', [0.4660, 0.6740, 0.1880, 0.5]);
% Experimental Kt
plot(sortedData.J, sortedData.Kt, 'o', 'MarkerFaceColor', [0, 0.4470, 0.7410], 'MarkerSize', 4, 'DisplayName', '$K_t$ (Experimental)', 'MarkerEdgeColor', [0, 0.4470, 0.7410]); 
% Polynomial fit for Kt
p_Kt = polyfit(sortedData.J, sortedData.Kt, 3);
Kt_poly = polyval(p_Kt, sortedData.J);
plot(sortedData.J, Kt_poly, '--', 'LineWidth', 1, 'DisplayName', '$K_t$ (Fit)', 'Color', [0, 0.4470, 0.7410]);

% Calculate R^2 for Kt fit
Kt_residuals = sortedData.Kt - Kt_poly;
Kt_SSresid = sum(Kt_residuals.^2);
Kt_SStotal = sum((sortedData.Kt - mean(sortedData.Kt)).^2);
Kt_R2 = 1 - (Kt_SSresid / Kt_SStotal);

% Print the equation and R^2 for Kt
fprintf('K_t Fit Equation: K_t = %.4f J^3 + %.4f J^2 + %.4f J + %.4f\n', p_Kt(1), p_Kt(2), p_Kt(3), p_Kt(4));
fprintf('K_t R^2 = %.4f\n\n', Kt_R2);

% Experimental Eta
plot(sortedData.J, sortedData.eta_exp, 'o', 'MarkerFaceColor', [0.4660, 0.6740, 0.1880], 'MarkerSize', 4, 'DisplayName', '$\eta$ (Experimental)', 'MarkerEdgeColor', [0.4660, 0.6740, 0.1880]);
% Polynomial fit for Efficiency
p_eta = polyfit(sortedData.J, sortedData.eta_exp, 3);
eta_poly = polyval(p_eta, sortedData.J);
% plot(sortedData.J, eta_poly, '--', 'LineWidth', 1, 'DisplayName', '$\eta$ (Fit)', 'Color', [0.4660, 0.6740, 0.1880]);

% Calculate R^2 for Eta fit
eta_residuals = sortedData.eta_exp - eta_poly;
eta_SSresid = sum(eta_residuals.^2);
eta_SStotal = sum((sortedData.eta_exp - mean(sortedData.eta_exp)).^2);
eta_R2 = 1 - (eta_SSresid / eta_SStotal);

% Print the equation and R^2 for Eta
fprintf('Eta Fit Equation: eta = %.2f J^3 + %.2f J^2 + %.2f J + %.2f\n', p_eta(1), p_eta(2), p_eta(3), p_eta(4));
fprintf('Eta R^2 = %.4f\n\n', eta_R2);

% Set Figure
ylabel('$K_t$, $\eta$', 'Interpreter', 'latex');
ylim([0, 1]);
ax.YColor = 'black';

% Right axis
yyaxis right;
% Reference Kq
plot(J_theo, K_Q_theo, '-', 'LineWidth', 2, 'DisplayName', '$K_q$ (Reference)', 'Color', [0.8500, 0.3250, 0.0980, 0.5]); 
hold on;
% Experimental Kq
plot(sortedData.J, sortedData.Kq, 'o', 'MarkerFaceColor', [0.8500, 0.3250, 0.0980], 'MarkerSize', 4, 'DisplayName', '$K_q$ (Experimental)', 'MarkerEdgeColor', [0.8500, 0.3250, 0.0980]); 

% Polynomial fit for Kq
p_Kq = polyfit(sortedData.J, sortedData.Kq, 3);
Kq_poly = polyval(p_Kq, sortedData.J);
plot(sortedData.J, Kq_poly, '--', 'LineWidth', 1, 'DisplayName', '$K_q$ (Fit)', 'Color', [0.8500, 0.3250, 0.0980]);

% Calculate R^2 for Kq fit
Kq_residuals = sortedData.Kq - Kq_poly;
Kq_SSresid = sum(Kq_residuals.^2);
Kq_SStotal = sum((sortedData.Kq - mean(sortedData.Kq)).^2);
Kq_R2 = 1 - (Kq_SSresid / Kq_SStotal);

% Print the equation and R^2 for Kq
fprintf('K_q Fit Equation: K_q = %.5f J^3 + %.5f J^2 + %.5f J + %.5f\n', p_Kq(1), p_Kq(2), p_Kq(3), p_Kq(4));
fprintf('K_q R^2 = %.4f\n\n', Kq_R2);

% Set right y-axis color to black
ax.YColor = 'black';
ylabel('$K_q$', 'Interpreter', 'latex');
ylim([0, 0.2]);

% Set up axis labels and title with LaTeX math notation
xlabel('$J$ (Advance Coefficient)', 'Interpreter', 'latex');
xlim([0, 0.773]);
title('Comparison of $K_t$, $K_q$, and $\eta$ vs $J$', 'Interpreter', 'latex');
legend('show', 'Interpreter', 'latex');
grid on;
