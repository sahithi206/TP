% probability.m
clear; clc;

% Data Used in Modelling 
load_means = [0.4526, 0.4526, 0.4526, 0.4715, 0.4715, 0.4715, 0.6552, 0.6552, 0.6552, 0.4919, 0.4919, 0.4919];
irradiance_data = [1.14, 2.3, 3.78, 5.08, 5.75, 6.27, 6.06, 5.25, 3.85, 2.51, 1.17, 0.9];

% Global variables
global gamma_Ct num_states num_solar_states;

%% Load modeling
mean_load = mean(load_means);
std_load = std(load_means);

num_states = 4; % Number of load states
load_min = max(mean_load - 3*std_load, 0); % Ensure non-negative
load_max = mean_load + 3*std_load;
state_bounds = linspace(load_min, load_max, num_states + 1);

% Compute Gamma Lt (load state probabilities)
gamma_Lt_values = zeros(num_states, 1);
fprintf('State Boundaries and Gamma Lt values:\n');
for i = 1:num_states
    state_start = state_bounds(i);
    state_end = state_bounds(i+1);
    gamma_Lt_values(i) = normcdf(state_end, mean_load, std_load) - normcdf(state_start, mean_load, std_load);
    fprintf('State %d: [%.4f, %.4f] -> Gamma Lt: %.4f\n', i, state_start, state_end, gamma_Lt_values(i));
end

% Replace zero/negative probabilities with tiny positive value
gamma_Lt_values(gamma_Lt_values <= 0) = eps;
gamma_Lt_values = gamma_Lt_values / sum(gamma_Lt_values); % Normalize

%% Solar irradiance modeling
min_val = min(irradiance_data);
max_val = max(irradiance_data);

% Handle case where all irradiance values are identical
if (max_val - min_val) == 0
    normalized_irradiance = zeros(size(irradiance_data));
else
    normalized_irradiance = (irradiance_data - min_val) / (max_val - min_val);
end

mean_normalized = mean(normalized_irradiance);
std_normalized = std(normalized_irradiance);

% Beta distribution parameters (ensure validity)
if std_normalized > 0 && (mean_normalized > 0) && (mean_normalized < 1)
    alpha = mean_normalized * ((mean_normalized*(1 - mean_normalized))/std_normalized^2 - 1);
    beta_param = (1 - mean_normalized) * ((mean_normalized*(1 - mean_normalized))/std_normalized^2 - 1);
else
    % Fallback to uniform distribution if parameters are invalid
    alpha = 1;
    beta_param = 1;
end

num_solar_states = 4; % Number of solar states
state_bounds_normalized = linspace(1, 10, num_solar_states + 1);
state_bounds_original = min_val + state_bounds_normalized * (max_val - min_val);
state_bounds_pairs = [state_bounds_original(1:end-1)', state_bounds_original(2:end)'];

% Compute Gamma St (solar state probabilities)
gamma_St_values = zeros(num_solar_states, 1);
for i = 1:num_solar_states
    lower = (state_bounds_pairs(i,1) - min_val) / (max_val - min_val);
    upper = (state_bounds_pairs(i,2) - min_val) / (max_val - min_val);
    gamma_St_values(i) = betacdf(upper, alpha, beta_param) - betacdf(lower, alpha, beta_param);
end

% Replace zero/negative probabilities with tiny positive value
gamma_St_values(gamma_St_values <= 0) = eps;
gamma_St_values = gamma_St_values / sum(gamma_St_values); % Normalize

%% Combined probabilities
gamma_Ct = gamma_Lt_values * gamma_St_values';
gamma_Ct = gamma_Ct / sum(gamma_Ct(:)); % Ensure normalization

% Display results
disp('Joint Probability Matrix (gamma_Ct):');
disp(gamma_Ct);