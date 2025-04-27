run('nlp.m');
run('powerLoss.m');
bus_numbers = 1:33;

global V_without_DG V_with_DG DG_buses x_opt fval ploss;

plotVoltageComparison(V_without_DG, V_with_DG, bus_numbers);
% Results
disp('V_with_DG');
disp(V_with_DG);
disp('V_without_DG');
disp(V_without_DG);
disp('Optimal DG sizes (MW):');
disp('Candidate Buses:');
fprintf('%d        ', DG_buses);
fprintf('\n');

disp('BM DGs:');
disp(x_opt(1:length(DG_buses))');
disp('PV systems:');
disp(x_opt(length(DG_buses)+1:end)');
disp(['Annual energy loss After Placing DG: ', num2str(fval), ' MWh']);
disp(['Annual energy loss Before Placing DG: ', num2str(ploss), ' MWh']);

function plotVoltageComparison(V_without_DG, V_with_DG, bus_numbers)
figure;
hold on;
grid on;

% Plot voltage profiles
plot(bus_numbers, V_without_DG(bus_numbers), '-', 'Color', [0 0.447 0.741], 'LineWidth', 2, 'DisplayName', 'Without DG');
plot(bus_numbers, V_with_DG(bus_numbers), '-', 'Color', [0.85 0.325 0.098], 'LineWidth', 2, 'DisplayName', 'With DGs');

% Add voltage limits
yline(0.95, 'k--', 'LineWidth', 1.2, 'DisplayName', 'Lower Limit (0.95 p.u.)');
yline(1.05, 'k--', 'LineWidth', 1.2, 'DisplayName', 'Upper Limit (1.05 p.u.)');

% Customize plot
xlabel('Bus Number', 'FontSize', 12);
ylabel('Voltage Magnitude (p.u.)', 'FontSize', 12);
title('Voltage Profile Comparison', 'FontSize', 14);

% Set axis limits and ticks
xlim([min(bus_numbers)-1 max(bus_numbers)+1]);
ylim([0.89 1.03]);
yticks(0.90:0.02:1.02);

% Add legend
legend('Location', 'best');
hold off;

% Save figure
saveas(gcf, 'voltage_comparison.png');
end

