clear; clc;
close all;
run('probability.m')
run('candidate.m')
global DG_buses voltage_limits S_B_min S_B_max I_max pf load_state_factors demand num_states mpc x0 num_solar_states gamma_Ct power_flow_results;
global V_without_DG V_with_DG x_opt fval;

% 1. Define system parameters
voltage_limits = [0.95, 1.05];     % PU voltage limits
S_B_min = 0.01;                     % Minimum DG size (MW)
S_B_max = 0.6;                        % Maximum DG size (MW)
I_max = 1; 
pf = 0.95;
load_state_factors = linspace(load_min, load_max, num_states); % Load factors
demand = sum(load_means);

results_without_DG = runpf(mpc);
V_without_DG = results_without_DG.bus(:, 8);  % Voltage magnitudes (p.u.)

% 3. Optimization setup
lb = repmat(S_B_min, 1, length(DG_buses)); % Lower bounds
ub = repmat(S_B_max, 1, length(DG_buses)); % Upper bounds
x0 = 0.305 * ones(length(DG_buses), 1);

% Pre-allocate storage for all power flow results
power_flow_results = struct();
power_flow_results.V = cell(num_states, num_solar_states);
power_flow_results.I = cell(num_states, num_solar_states);
power_flow_results.P_flow = cell(num_states, num_solar_states);
power_flow_results.Q_flow = cell(num_states, num_solar_states);

function [V, I, P_flow, Q_flow] = SolvePowerFlow(x, load_state, solar_state)
    global load_state_factors DG_buses mpc V_with_DG power_flow_results;
    
    % Check if we already have these results stored
    if ~isempty(power_flow_results.V{load_state, solar_state})
        V = power_flow_results.V{load_state, solar_state};
        I = power_flow_results.I{load_state, solar_state};
        P_flow = power_flow_results.P_flow{load_state, solar_state};
        Q_flow = power_flow_results.Q_flow{load_state, solar_state};
        return;
    end
    
    mpc_mod = mpc;
    load_factor = load_state_factors(load_state);
    mpc_mod.bus(:, 3:4) = mpc.bus(:, 3:4) * load_factor;

    % Retain only the slack generator initially
    mpc_mod.gen = mpc_mod.gen(1,:);
    mpc_mod.gencost = mpc_mod.gencost(1,:);

    % Vectorized DG addition
    newGens = zeros(length(DG_buses), 21);
    newGens(:, 1) = DG_buses;
    newGens(:, 2) = x;
    newGens(:, 6) = 1.0;
    newGens(:, 7) = 100;
    newGens(:, 8) = 1;
    newGens(:, 9) = x;
    newGens(:, [4,5]) = repmat([100, -100], length(DG_buses), 1);
    
    mpc_mod.gen = [mpc_mod.gen; newGens];
    mpc_mod.gencost = repmat([2, 0, 0, 3, 0.1, 5, 0], length(DG_buses)+1, 1);

    results = runpf(mpc_mod);
    V_with_DG = results.bus(:, 8);  
    V = results.bus(:, 8);
    V_from = results.bus(results.branch(:,1), 8);
    V_to = results.bus(results.branch(:,2), 8);
    Z = results.branch(:,3) + 1j*results.branch(:,4);
    I = abs((V_from - V_to) ./ Z);
    P_flow = results.branch(:, 14);
    Q_flow = results.branch(:, 15);
    
    % Store results for future use
    power_flow_results.V{load_state, solar_state} = V;
    power_flow_results.I{load_state, solar_state} = I;
    power_flow_results.P_flow{load_state, solar_state} = P_flow;
    power_flow_results.Q_flow{load_state, solar_state} = Q_flow;
end

function P_Loss = MinPLoss(x)
    global num_states num_solar_states mpc gamma_Ct power_flow_results V_sum;
    P_Loss = 0;
    V_sum=0;
    for i = 1:num_states
        for j = 1:num_solar_states
            [V, I] = SolvePowerFlow(x, i, j);
            Rij = mpc.branch(:, 3);
            P_Loss_ij = sum(I.^2 .* Rij) * mpc.baseMVA;
            P_Loss = P_Loss + P_Loss_ij * gamma_Ct(i, j) * 8760;
        end
    end
end

function [c, ceq] = PowerConstraintsMulti(x)
    global num_states num_solar_states mpc gamma_Ct power_flow_results voltage_limits I_max;

    % Clear previous power flow results when x changes
    power_flow_results.V = cell(num_states, num_solar_states);
    power_flow_results.I = cell(num_states, num_solar_states);
    power_flow_results.P_flow = cell(num_states, num_solar_states);
    power_flow_results.Q_flow = cell(num_states, num_solar_states);

    % Pre-allocate constraints
    c = [];
    ceq = [];

    for i = 1:num_states
        for j = 1:num_solar_states
            [V, I] = SolvePowerFlow(x, i, j);

            % Voltage constraints
            voltage_upper = V - voltage_limits(2);
            voltage_lower = voltage_limits(1) - V;

            % Current constraints
            current_limits = abs(I) - I_max;

            % Load calculation
            total_load = sum(mpc.bus(:,3));
            total_gen = sum(x);
            gen_demand_constraint = total_gen - total_load;

            % Combine constraints
            constraints = [voltage_upper; voltage_lower; current_limits; gen_demand_constraint];
            c = [c; constraints];
        end
    end
end


options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'interior-point', ...
                     'MaxFunctionEvaluations', 5000, 'UseParallel', true);
minploss_fun = @(x) MinPLoss(x);
nonlcon_fun = @(x) PowerConstraintsMulti(x);
[x_opt, fval] = fmincon(minploss_fun, x0, [], [], [], [], lb, ub, nonlcon_fun, options);

disp('Optimal DG sizes (MW):');
disp(x_opt');
disp(['Annual energy loss: ', num2str(fval), ' MWh']);

total_voltage_with_DG = 0;
num_entries_with_DG = 0;

for i = 1:num_states
    for j = 1:num_solar_states
        [V, ~, ~, ~] = SolvePowerFlow(x_opt, i, j); 
        total_voltage_with_DG = total_voltage_with_DG + V;
        num_entries_with_DG = num_entries_with_DG + 1;
    end
end

V_with_DG = total_voltage_with_DG / num_entries_with_DG