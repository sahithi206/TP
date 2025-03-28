%%load and irradiance states
% Number of load states and irradiance states
num_load_states = 10;
num_irrad_states = 9;

% Load and irradiance data for each month
load_data = [0.4526, 0.4526, 0.4526, 0.4715, 0.4715, 0.4715, ...
             0.6552, 0.6552, 0.6552, 0.4919, 0.4919, 0.4919];
irradiance_data = [1.14, 2.3, 3.78, 5.08, 5.75, 6.27, ...
                   6.06, 5.25, 3.85, 2.51, 1.17, 0.9];

% Normalize load and irradiance data
norm_load = rescale(load_data, 0, 1);
norm_irradiance = rescale(irradiance_data, 0, 1);

% Generate discrete states using clustering (K-means)
[idx_load, load_centroids] = kmeans(norm_load', num_load_states);
[idx_irradiance, irradiance_centroids] = kmeans(norm_irradiance', num_irrad_states);

% Assign load states and irradiance states to months
load_states = load_centroids(idx_load);
irradiance_states = irradiance_centroids(idx_irradiance);

% Assume an example bus system with N buses
num_buses = 5;  % Modify as per the system
bus_load_distribution = zeros(num_buses, num_load_states);

% Distribute load among buses proportionally
for i = 1:num_load_states
    total_load = load_states(i);
    bus_load_distribution(:, i) = total_load * rand(num_buses, 1);
    bus_load_distribution(:, i) = bus_load_distribution(:, i) / sum(bus_load_distribution(:, i)) * total_load;
end

% Display Results
disp('Load states:');
disp(load_states');
disp('Irradiance states:');
disp(irradiance_states');
disp('Bus Load Distribution:');
disp(bus_load_distribution);
