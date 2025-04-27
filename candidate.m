global mpc;

mpc = loadcase('case33bw');
results=runpf(mpc);
V_base = results.bus(:, 8);  
mpc.bus(:, 3) = mpc.bus(:, 3) * 1.072;  % Increase active power demand (P)
mpc.bus(:, 4) = mpc.bus(:, 4) * 1.072;  % Increase reactive power demand (Q)

% Solve Power Flow for Load Growth Case
results_after = runpf(mpc);
V_after = results_after.bus(:, 8);

% Compute Voltage Regulation Index (VRI)
VRI = (V_base - V_after)./ V_after;
disp("VRI"+VRI);

threshold= max(VRI)*0.95;
candidate_buses = find(VRI >threshold);

DG_buses=candidate_buses;  
candidates_matrix = [candidate_buses, VRI(candidate_buses)];

% Display Results
disp("Candidate Buses :");
disp(candidates_matrix);