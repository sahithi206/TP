clear; clc;
close all;
run('probability.m')
disp("LOad States");
disp(num_states)
global DG_buses voltage_limits S_B_min S_B_max I_max   pf  load_state_factors;
% 1. Define system parameters
DG_buses = [16, 17, 18, 32, 33];  % DG locations
voltage_limits = [0.95, 1.05];     % PU voltage limits
S_B_min = 0.1;                     % Minimum DG size (MW)
S_B_max = 2.0;                     % Maximum DG size (MW)
I_max = 300;                       % Thermal limit (A)
pf = 0.95;
load_state_factors = linspace(load_min, load_max, num_states); % Load factors

mpc = loadcase('case33bw');
disp('Generator data (mpc.gen):');
disp(mpc.gen);
disp('Generator cost data (mpc.gencost):');
disp(mpc.gencost);

% 3. Optimization setup
lb = repmat(S_B_min, 1, length(DG_buses)); % Lower bounds
ub = repmat(S_B_max, 1, length(DG_buses)); % Upper bounds
x0=[0.5,0.5,0.5,0.5,0.5]
function [V, I, P_flow, Q_flow] = SolvePowerFlow(x, load_state, solar_state)
    global load_state_factors DG_buses mpc;

    % Copy the base case data
    mpc_mod = mpc;
    
    % Adjust load levels
    load_factor = load_state_factors(load_state);
    mpc_mod.bus(:, 3) = mpc.bus(:, 3) * load_factor; % P
    mpc_mod.bus(:, 4) = mpc.bus(:, 4) * load_factor; % Q

    % Clear existing generators (keep only slack bus)
   mpc_mod.gen = mpc_mod.gen(1,:);
    mpc_mod.gencost = mpc_mod.gencost(1,:);

    % Add DGs as generators at specified buses
    for i = 1:length(DG_buses)
        if DG_buses(i) == 0
            continue;
        end
        
        newGen = zeros(1, 21);
        newGen(1, 1) = DG_buses(i);       % bus number
        newGen(1, 2) = x(i);              % PG
        newGen(1, 3) = 0;                 % QG
        newGen(1, 4) = 100;               % QMAX
        newGen(1, 5) = -100;              % QMIN
        newGen(1, 6) = 1.0;              % VG
        newGen(1, 7) = 100;               % MBASE
        newGen(1, 8) = 1;                 % GEN STATUS
        newGen(1, 9) = x(i);              % PMAX
        newGen(1, 10) = 0;                % PMIN
        % Remaining columns can stay 0

        % Create corresponding gencost entry (polynomial cost)
          if isempty(mpc_mod.gen)
            mpc_mod.gen = newGen; % Assign directly if empty
          elseif DG_buses(i)==0
               break
           else
             mpc_mod.gen = [mpc_mod.gen; newGen];

            end
    end
    % Polynomial cost (model=2) with 3 coefficients (n=3)
mpc_mod.gencost = [
    2, 0, 0, 3, 0.1, 5, 0;  % Generator 1 cost: 0.1 + 5*P + 0*P²
    2, 0, 0, 3, 0.1, 5, 0;  % Generator 2 cost
    2, 0, 0, 3, 0.1, 5, 0;  % Generator 3 cost
    2, 0, 0, 3, 0.1, 5, 0;  % Generator 4 cost
    2, 0, 0, 3, 0.1, 5, 0;  % Generator 5 cost
    2, 0, 0, 3, 0.1, 5, 0;  % Generator 6 cost
];
    disp(mpc_mod);
    % Run Power Flow
    results = runpf(mpc_mod);
    % Extract results
    V = results.bus(:, 8);
    I = abs(results.branch(:, 14) + 1j * results.branch(:, 15)) ./ results.branch(:, 3);
    P_flow = results.branch(:, 14);
    Q_flow = results.branch(:, 15);
end

function P_Loss = MinPLoss(x)
global num_states num_solar_states mpc gamma_Ct;
% Fixed DG locations
DG_buses = [16, 17, 18, 32, 33];

P_Loss = 0; % Initialize total loss

% Loop over all load and solar states
disp("States:");
disp(num_states);
for i = 1:num_states
    for j = 1:num_solar_states
        % Solve power flow for (load state i, solar state j)
        [V, I,P_flow,Q_flow] = SolvePowerFlow(x, i, j); % Updated function

        % Extract resistance values
        Rij = mpc.branch(:, 3);  % Assuming column 3 is resistance


        % Compute power loss for this state
        P_Loss_ij = sum(sum((I.^2) .* Rij));

        % Weight by joint probability γ(Ct) for (i, j)
        P_Loss = P_Loss + P_Loss_ij * gamma_Ct(i, j) * 8760;
    end
end
end

function [c, ceq] = PowerFlowConstraints(x, load_state, solar_state)
    global pf voltage_limits I_max
    
    % Solve power flow with given DG outputs x
    [V, I, P_flow, Q_flow] = SolvePowerFlow(x, load_state, solar_state);
    
    % Inequality constraints (c ≤ 0):
    c = [
        V - voltage_limits(2);      % Ensure V ≤ Vmax (for all buses)
        voltage_limits(1) - V;      % Ensure V ≥ Vmin (for all buses)
        I - I_max;                  % Ensure I ≤ Imax (for all branches)
    ];
    
    % Note: x constraints (S_B_min ≤ x ≤ S_B_max) should be handled via lb/ub in fmincon
    % rather than here in nonlinear constraints for better numerical performance
    
    % Equality constraints (ceq = 0):
    ceq = [];  % No equality constraints
end

options = optimoptions('fmincon',...
    'Display', 'iter',...
    'Algorithm', 'interior-point',...
    'MaxFunctionEvaluations', 5000);

% 4. Create function handles
minploss_fun = @(x) MinPLoss(x);
nonlcon_fun = @(x) PowerFlowConstraints(x);

% 5. Run optimization
[x_opt, fval] = fmincon(minploss_fun, x0, [], [], [], [], lb, ub, nonlcon_fun, options);

% 6. Display results
disp('Optimal DG sizes (MW):');
disp(x_opt');
disp(['Annual energy loss: ', num2str(fval), ' MWh']);