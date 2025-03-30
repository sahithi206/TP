clear; clc;
close all;
run('probability.m')
global DG_buses voltage_limits S_B_min S_B_max I_max pf load_state_factors demand num_states num_solar_states mpc gamma_Ct;
global V_without_DG  V_with_DG;
% 1. Define system parameters
DG_buses = [16, 17, 18, 32, 33,0];  % DG locations
voltage_limits = [0.95, 1.05];     % PU voltage limits
S_B_min = 0.01;                     % Minimum DG size (MW)
S_B_max = 2;                     % Maximum DG size (MW)
pf = 0.95;
load_state_factors = linspace(load_min,load_max, num_states); % Load factors
demand = sum(load_means);
mpc = loadcase('case33bw');

results_without_DG = runpf(mpc);
V_without_DG = results_without_DG.bus(:, 8);  % Voltage magnitudes (p.u.)

% 3. Optimization setup
lb = repmat(S_B_min, 1, length(DG_buses)); % Lower bounds
ub = repmat(S_B_max, 1, length(DG_buses)); % Upper bounds
x0 = [0.04,0.04,0.04,0.04,0.04]; % Initial guess within bounds

function [V, I, P_flow, Q_flow] = SolvePowerFlow(x, load_state, solar_state)
    global load_state_factors DG_buses mpc V_with_DG;
    mpc_mod = mpc;
    load_factor = load_state_factors(load_state);
    mpc_mod.bus(:, 3) = mpc.bus(:, 3) * load_factor;
    mpc_mod.bus(:, 4) = mpc.bus(:, 4) * load_factor;
    
    mpc_mod.gen = mpc_mod.gen(1,:);
    mpc_mod.gencost = mpc_mod.gencost(1,:);
    
    for i = 1:length(DG_buses)
        if DG_buses(i) == 0
            continue;
        end
        
        newGen = zeros(1, 21);
        newGen(1, 1) = DG_buses(i);
        newGen(1, 2) = x(i);
        newGen(1, 3) = 0;
        newGen(1, 4) = 100;
        newGen(1, 5) = -100;
        newGen(1, 6) = 1.0;
        newGen(1, 7) = 100;
        newGen(1, 8) = 1;
        newGen(1, 9) = x(i);
        newGen(1, 10) = 0;
        
        if isempty(mpc_mod.gen)
            mpc_mod.gen = newGen;
        elseif DG_buses(i) == 0
            break;
        else
            mpc_mod.gen = [mpc_mod.gen; newGen];
        end
    end
    
    mpc_mod.gencost = repmat([2, 0, 0, 3, 0.1, 5, 0], length(DG_buses), 1);
    results = runpf(mpc_mod);
   V_with_DG = results.bus(:, 8);  % Voltage magnitudes (p.u.)
    V = results.bus(:, 8);
    I = abs(results.branch(:, 14) + 1j * results.branch(:, 15)) ./ results.branch(:, 3);
    P_flow = results.branch(:, 14);
    Q_flow = results.branch(:, 15);
end

function P_Loss = MinPLoss(x)
    global num_states num_solar_states mpc gamma_Ct;
    DG_buses = [16, 17, 18, 32, 33];
    P_Loss = 0;
    for i = 1:num_states
        for j = 1:num_solar_states
            [V, I, P_flow, Q_flow] = SolvePowerFlow(x, i, j);
            Rij = mpc.branch(:, 3);
            P_Loss_ij = sum(sum((I.^2) .* Rij));
            P_Loss = P_Loss + P_Loss_ij * gamma_Ct(i, j) * 8760;
        end
    end
end
function [c, ceq] = PowerConstraintsMulti(x)
    global num_states num_solar_states mpc gamma_Ct;
    
    % Initialize constraints
    c = []; % Empty array to store inequality constraints
    ceq = []; % No equality constraints
    
    % Loop over all load and solar states
    for i = 1:num_states
        for j = 1:num_solar_states
            % Get constraints for this (load, solar) state
            [constraints, constraintseq] = PowerFlowConstraints(x, i, j);
            
            % Concatenate inequality constraints
            c = [c; constraints];
            
            % Concatenate equality constraints (if needed)
            ceq = [ceq; constraintseq];
        end
    end
end


function [constraints, constraintseq] = PowerFlowConstraints(x, load_state, solar_state)
    global pf voltage_limits I_max demand;
    [V, I, P_flow, Q_flow] = SolvePowerFlow(x, load_state, solar_state);
    constraints = [V - voltage_limits(2); voltage_limits(1) - V];
    constraintseq = [];
end

options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'interior-point', 'MaxFunctionEvaluations', 5000);
minploss_fun = @(x) MinPLoss(x);
nonlcon_fun = @(x) PowerConstraintsMulti(x);
[x_opt, fval] = fmincon(minploss_fun, x0, [], [], [], [], lb, ub, nonlcon_fun, options);

disp('Optimal DG sizes (MW):');
disp(x_opt');
disp(['Annual energy loss: ', num2str(fval), ' MWh']);


function plotVoltageComparison(V_without_DG, V_with_DG, bus_numbers)
    % Create figure
    figure;
    hold on;
    grid on;
    
    % Plot voltage profiles
   % Plot voltage profiles with smooth lines and no markers
plot(bus_numbers, V_without_DG(bus_numbers), '-', 'Color', [0 0.447 0.741], ...  % Blue
    'LineWidth', 2, 'DisplayName', 'Without DG');
plot(bus_numbers, V_with_DG(bus_numbers), '-', 'Color', [0.85 0.325 0.098], ...   % Red
    'LineWidth', 2, 'DisplayName', 'With DGs');
    
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
% Select buses to plot (e.g., all buses or specific ones)
bus_numbers = 1:33;  % All buses
% Or for odd buses: bus_numbers = 3:2:33;

% Generate the plot
plotVoltageComparison(V_without_DG, V_with_DG, bus_numbers);