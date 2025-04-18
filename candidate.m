mpc = loadcase('case33bw');
results=runpf(mpc);

V_base = results.bus(:, 8);  
mpc.bus(:, 3) = mpc.bus(:, 3) * 1.072;  % P incr
mpc.bus(:, 4) = mpc.bus(:, 4) * 1.072;  % Q incr

results_after = runpf(mpc);
V_after = results_after.bus(:, 8);

VRI = (V_base - V_after)./ V_after;
disp("VRI"+VRI);
candidate_buses = find(VRI >0.00718);

DG_buses=candidate_buses;  
disp("Buses");
disp(candidate_buses)
candidates_matrix = [candidate_buses, VRI(candidate_buses)];

disp("Candidate Buses Where VRI > 0.0072:");
disp(candidates_matrix);