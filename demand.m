% Define load and irradiance data
load_data = [0.4526, 0.4526, 0.4526, 0.4715, 0.4715, 0.4715, ...
             0.6552, 0.6552, 0.6552, 0.4919, 0.4919, 0.4919];
irradiance_data = [1.14, 2.3, 3.78, 5.08, 5.75, 6.27, ...
                   6.06, 5.25, 3.85, 2.51, 1.17, 0.9];

% Define number of states
num_load_states = 10;
num_irrad_states = 9;s

% Compute actual min and max of load data
load_min = floor(min(load_data)); % Round down
load_max = ceil(max(load_data));  % Round up

% Define bin size
interval_width = (load_max - load_min) / num_load_states;

% Initialize cell array to store intervals
load_states = cell(1, num_load_states);
% Generate Load Intervals
for i = 1:num_load_states
    lower_bound = load_min + (i - 1) * interval_width;
    upper_bound = lower_bound + interval_width;

    % Assign interval to the cell array
    load_states{i} = [lower_bound, upper_bound];

    % Condition to avoid exceeding load_max
    if upper_bound >= load_max
        break;
    end
end

% Display results
disp('Load states (MW):');
disp(load_states');
