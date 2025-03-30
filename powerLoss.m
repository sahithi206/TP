clear; clc;
close all;
run('probability.m')
global DG_buses voltage_limits S_B_min S_B_max I_max   pf  load_state_factors demand;
% 1. Define system parameters
voltage_limits = [0.95, 1.05];     % PU voltage limits
S_B_min = 0.1;                     % Minimum DG size (MW)
S_B_max = 2.0;                     % Maximum DG size (MW)
I_max = 300;                       % Thermal limit (A)
pf = 0.95;
load_state_factors = linspace(load_min, load_max, num_states); % Load factors
demand=sum(load_means);
mpc = loadcase('case33bw');

function [V, I, P_flow, Q_flow] = SolvePowerFlow(load_state, solar_state)
    global load_state_factors DG_buses mpc;

    % Copy the base case data
    mpc_mod = mpc;
    
    % Adjust load levels
    load_factor = load_state_factors(load_state);
    mpc_mod.bus(:, 3) = mpc.bus(:, 3) * load_factor; % P
    mpc_mod.bus(:, 4) = mpc.bus(:, 4) * load_factor; % Q

    % Clear existing generators (keep only slack bus)
   mpc_mod.gen = mpc_mod.gen(1,:);
    mpc_mod.gencost = mpc_mod.gencost(1,:);

    % Run Power Flow
    results = runpf(mpc_mod);
    % Extract results
    V = results.bus(:, 8);
    I = abs(results.branch(:, 14) + 1j * results.branch(:, 15)) ./ results.branch(:, 3);
    P_flow = results.branch(:, 14);
    Q_flow = results.branch(:, 15);
end

function P_Loss = MinPLoss()
    global num_states num_solar_states mpc gamma_Ct;
    
    P_Loss = 0; % Initialize total loss
    
    % Loop over all load and solar states
    for i = 1:num_states
        for j = 1:num_solar_states
            % Solve power flow for (load state i, solar state j)
            [V, I, P_flow, Q_flow] = SolvePowerFlow(i, j); 
            
            % Extract resistance values
            Rij = mpc.branch(:, 3);  % Assuming column 3 contains resistance
            
            % Compute power loss for this state
            P_Loss_ij = sum((I.^2) .* Rij);  
            
            % Weight by joint probability γ(Ct) for (i, j)
            P_Loss = P_Loss + P_Loss_ij * gamma_Ct(i, j) * 8760;
        end
    end
end

P_Loss_no_DG = MinPLoss();
disp(['Annual energy loss without DGs: ', num2str(P_Loss_no_DG)]);