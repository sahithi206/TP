run('demand.m'); 
V_base = results.bus(:, 8);  

mpc.bus(:, 3) = mpc.bus(:, 3) * 1.072;  % Increase active power demand (P)
mpc.bus(:, 4) = mpc.bus(:, 4) * 1.072;  % Increase reactive power demand (Q)

% Solve Power Flow for Load Growth Case
results_after = runpf(mpc);
V_after = results_after.bus(:, 8);  % Extract updated voltage magnitudes

% Compute Voltage Regulation Index (VRI)
VRI = (V_base - V_after) ./ V_after;

disp("VRI"+VRI);

candidate_buses = find(VRI > 0.000798);  % Adjusted threshold
disp("Buses");
disp(candidate_buses)

% Store Selected Buses in a Matrix
candidates_matrix = [candidate_buses, VRI(candidate_buses)];

% Display Results
disp("Candidate Buses Where VRI > 0.000798:");
disp(candidates_matrix);
