define_constants;
mpc = loadcase('case33bw');  % Load IEEE 33-bus case

% Given: load_data (12 months of load demand)
annual_avg_load = mean(load_data); % Compute annual average load

% Extract original IEEE-33 bus load values
original_PD = mpc.bus(:, PD); 
original_QD = mpc.bus(:, QD);

% Compute total original IEEE-33 load demand
total_original_PD = sum(original_PD);

% Compute scaling factor based on new annual average load
scaling_factor = annual_avg_load / total_original_PD;

% Scale original IEEE-33 loads while maintaining proportions
P_load = original_PD * scaling_factor; 
Q_load = original_QD * scaling_factor;

% Assign to MATPOWER case struct
mpc.bus(:, PD) = P_load;
mpc.bus(:, QD) = Q_load;

% Run Power Flow Analysis
results = runpf(mpc);

% Display Results
disp('Power Flow Analysis Results:');

