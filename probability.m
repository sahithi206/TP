run('load_states.m'); 
disp(load_states);

mean_load = mean(cellfun(@(x) mean(x), load_states));
std_load = std(cellfun(@(x) mean(x), load_states));

% Compute Gamma Lt for each load state
gamma_Lt_values = zeros(num_load_states, 1);
for i = 1:num_load_states
    state_start = load_states{i}(1);
    state_end = load_states{i}(2);
    gamma_Lt_values(i) = normcdf(state_end, mean_load, std_load) - normcdf(state_start, mean_load, std_load);
end

normalized_irradiance = (irradiance_data - irrad_min) / (irrad_max - irrad_min);
mean_normalized = mean(normalized_irradiance);
std_normalized = std(normalized_irradiance);

% Compute Beta distribution parameters
if std_normalized > 0
    alpha = mean_normalized * ((mean_normalized * (1 - mean_normalized)) / std_normalized^2 - 1);
    beta_param = (1 - mean_normalized) * ((mean_normalized * (1 - mean_normalized)) / std_normalized^2 - 1);
else
    alpha = 1; 
    beta_param = 1; % Default values to prevent division by zero
end

% Compute Gamma St values (probabilities for each irradiance state)
num_solar_states = 10;
state_bounds_normalized = linspace(0, 1, num_solar_states + 1);
state_bounds_original = min_val + state_bounds_normalized * (max_val - min_val);
state_bounds_pairs = [state_bounds_original(1:end-1)', state_bounds_original(2:end)'];

% Compute Gamma St
gamma_St_values = zeros(size(state_bounds_pairs, 1), 1);
for i = 1:size(state_bounds_pairs, 1)
    lower_bound = (state_bounds_pairs(i,1) - min_val) / (max_val - min_val);
    upper_bound = (state_bounds_pairs(i,2) - min_val) / (max_val - min_val);
    gamma_St_values(i) = betacdf(upper_bound, alpha, beta_param) - betacdf(lower_bound, alpha, beta_param);
end

disp(gamma_St_values);
% Compute joint probability
gamma_Ct = gamma_Lt_values * gamma_St_values';

disp(gamma_Ct)