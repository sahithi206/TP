load_means = [0.4526, 0.4526, 0.4526, ...  
              0.4715, 0.4715, 0.4715, ...  
              0.6552, 0.6552, 0.6552, ...  
              0.4919, 0.4919, 0.4919];     
irradiance_data = [1.14, 2.3, 3.78, 5.08, 5.75, 6.27, ...
                   6.06, 5.25, 3.85, 2.51, 1.17, 0.9];
global gamma_Ct  num_states num_solar_states;
% Compute combined mean and std deviation
mean_load = mean(load_means);
std_load = std(load_means);

% Define 9 load states (equally spaced)
num_states = 2;
load_min = mean_load - 3 * std_load;
load_max = mean_load + 3 * std_load;
state_bounds = linspace(load_min, load_max, num_states + 1);

% Compute Gamma Lt for each state
gamma_Lt_values = zeros(num_states, 1);

fprintf('State Boundaries and Gamma Lt values:\n');
for i = 1:num_states
    state_start = state_bounds(i);
    state_end = state_bounds(i+1);
    
    % Corrected formula for Gamma Lt
    gamma_Lt_values(i) = normcdf(state_end, mean_load, std_load) - normcdf(state_start, mean_load, std_load);
    
    fprintf('State %d: [%.4f, %.4f] -> Gamma Lt: %.4f\n', ...
            i, state_start, state_end, gamma_Lt_values(i));
end
irrad_min=min(irradiance_data);
irrad_max=max(irradiance_data);
normalized_irradiance = (irradiance_data - irrad_min) / (irrad_max - irrad_min);
mean_normalized = mean(normalized_irradiance);
std_normalized = std(normalized_irradiance);

% Compute Beta distribution parameters
if std_normalized > 0
    alpha = mean_normalized * ((mean_normalized * (1 - mean_normalized)) / std_normalized^2 - 1);
    beta_param = (1 - mean_normalized) * ((mean_normalized * (1 - mean_normalized)) / std_normalized^2 - 1);
else
    alpha = 1; 
    beta_param = 1; 
end
min_val=irrad_min;
max_val=irrad_max;

num_solar_states = 2;
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
