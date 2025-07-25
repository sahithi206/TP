clear; clc;
close all;
run('probability.m')
run('candidate.m')
global  load_state_factors demand num_states mpc num_solar_states  gamma_Ct ploss ;

load_state_factors = linspace(load_min,load_max, num_states);

function [V, I, P_flow, Q_flow] = SolvePowerFlow( load_state,solar_state)
global load_state_factors  mpc ;
mpc_mod = mpc;
load_factor = load_state_factors(load_state);
mpc_mod.bus(:, 3) = mpc.bus(:, 3) * load_factor;
mpc_mod.bus(:, 4) = mpc.bus(:, 4) * load_factor;
results = runpf(mpc_mod);
V = results.bus(:, 8);
V_from = results.bus(results.branch(:,1), 8);
V_to = results.bus(results.branch(:,2), 8);
Z = results.branch(:,3) + 1j*results.branch(:,4);
I = abs((V_from - V_to) ./ Z);
P_flow = results.branch(:, 14);
Q_flow = results.branch(:, 15);
end

function P_Loss = MinPLoss()
global num_states num_solar_states mpc gamma_Ct  ;
P_Loss = 0;
for i = 1:num_states
    for j = 1:num_solar_states
        [V, I, P_flow, Q_flow] = SolvePowerFlow(i, j);
        Rij = mpc.branch(:, 3);
        P_Loss_ij = sum(I.^2 .* Rij) * mpc.baseMVA;
        P_Loss = P_Loss + P_Loss_ij * gamma_Ct(i, j) * 8760;
    end
end
end
ploss=MinPLoss();

disp(ploss);

