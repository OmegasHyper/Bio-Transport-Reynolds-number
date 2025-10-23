%% Reynolds Number and Flow Regime Analysis
% SBEG201 Biotransport - Assignment 01
% Complete analysis of Reynolds number effects on fluid flow

clear; clc; close all;

%% Part 1: Reynolds Number Calculator & Regime Classification

fprintf('=== PART 1: Reynolds Number Calculator & Regime Classification ===\n\n');

% Define fluid properties database
fluids = {
    % Name, Density (kg/m³), Dynamic Viscosity (Pa·s), Type
    'Water at 5°C', 1000, 1.518e-3, 'Newtonian';
    'Water at 20°C', 998, 1.002e-3, 'Newtonian';
    'Water at 40°C', 992, 0.653e-3, 'Newtonian';
    'Air at 20°C', 1.204, 1.825e-5, 'Newtonian';
    'Blood (37°C)', 1060, 3.5e-3, 'Non-Newtonian';
    'CSF (37°C)', 1007, 0.8e-3, 'Newtonian';
    'Glycerin (20°C)', 1260, 1.410, 'Newtonian';
    'Honey (20°C)', 1420, 10.0, 'Newtonian'
};

% Pipe configurations
pipe_diameters = [0.005, 0.01, 0.025, 0.05, 0.1]; % meters
velocities = [0.01, 0.1, 0.5, 1.0, 2.0, 5.0]; % m/s

% Function to calculate Reynolds number
calculate_Re = @(rho, V, D, mu) (rho * V * D) / mu;

% Function to classify flow regime
classify_regime = @(Re) ...
    categorical(Re < 2000, [true, false], {'Laminar', 'Turbulent'});

% Generate and display results
results_table = table();
counter = 1;

fprintf('Calculating Reynolds numbers for various configurations...\n\n');

for fluid_idx = 1:size(fluids, 1)
    fluid_name = fluids{fluid_idx, 1};
    rho = fluids{fluid_idx, 2};
    mu = fluids{fluid_idx, 3};
    fluid_type = fluids{fluid_idx, 4};
    
    for diam_idx = 1:length(pipe_diameters)
        D = pipe_diameters(diam_idx);
        
        for vel_idx = 1:length(velocities)
            V = velocities(vel_idx);
            
            % Calculate Reynolds number
            Re = calculate_Re(rho, V, D, mu);
            
            % Classify regime
            if Re < 2000
                regime = 'Laminar';
            elseif Re < 4000
                regime = 'Transitional';
            else
                regime = 'Turbulent';
            end
            
            % Store results
            results_table.Fluid(counter) = {fluid_name};
            results_table.Density(counter) = rho;
            results_table.Viscosity(counter) = mu;
            results_table.Diameter(counter) = D;
            results_table.Velocity(counter) = V;
            results_table.Reynolds(counter) = Re;
            results_table.Regime(counter) = {regime};
            results_table.Type(counter) = {fluid_type};
            
            counter = counter + 1;
        end
    end
end

% -----------------------------------------------------------------
% --- Modified display logic to ensure diversity in fluids and systems ---
% -----------------------------------------------------------------
fprintf('Sample Results (Selected to show fluid and regime variation):\n');
selected_indices = [];
unique_fluids = unique(results_table.Fluid);

% 1. Ensure the first entry of each fluid appears (usually Laminar)
for i = 1:length(unique_fluids)
    fluid_name = unique_fluids{i};
    match_index = find(strcmp(results_table.Fluid, fluid_name), 1, 'first');
    if ~isempty(match_index)
        selected_indices = [selected_indices; match_index];
    end
end

% 2. Select additional entries for key fluids (Blood, Water, Air) to show regime transition
fluids_for_transition = {'Water at 20°C', 'Blood (37°C)', 'Air at 20°C'};

for i = 1:length(fluids_for_transition)
    current_fluid = fluids_for_transition{i};
    
    % A. Search for the Transitional case
    trans_index = find(strcmp(results_table.Fluid, current_fluid) & ...
                       strcmp(results_table.Regime, 'Transitional'), 1, 'first');
    if ~isempty(trans_index) && ~ismember(trans_index, selected_indices)
        selected_indices = [selected_indices; trans_index];
    end

    % B. Search for the Turbulent case
    turb_index = find(strcmp(results_table.Fluid, current_fluid) & ...
                      strcmp(results_table.Regime, 'Turbulent'), 1, 'first');
    if ~isempty(turb_index) && ~ismember(turb_index, selected_indices)
        selected_indices = [selected_indices; turb_index];
    end
end

% 3. Sort the selected rows alphabetically by fluid name
selected_data = results_table(selected_indices, :);
selected_data = sortrows(selected_data, 'Fluid'); % Sort fluids alphabetically

% Display the final table
if ~isempty(selected_data)
    disp(selected_data);
else
    fprintf('No results to display.\n');
end

% Summary statistics
fprintf('\nFlow Regime Summary:\n');
regime_summary = groupsummary(results_table, 'Regime');
disp(regime_summary);


%% Part 2: Flow Profile Visualization

fprintf('\n=== PART 2: Flow Profile Visualization ===\n\n');

% Select representative cases for visualization
cases_to_plot = [
    50, 2000;    % Laminar case
    150, 5000;   % Turbulent case
    80, 3500     % Transitional case
];

figure('Position', [100, 100, 1200, 800]);
sgtitle('Velocity Profiles for Different Flow Regimes', 'FontSize', 14, 'FontWeight', 'bold');

for case_idx = 1:size(cases_to_plot, 1)
    row_idx = cases_to_plot(case_idx, 1);
    Re_target = cases_to_plot(case_idx, 2);
    
    % Find case closest to target Re
    [~, closest_idx] = min(abs(results_table.Reynolds - Re_target));
    case_data = results_table(closest_idx, :);
    
    % Pipe geometry
    R = case_data.Diameter / 2; % Pipe radius
    r = linspace(-R, R, 100); % Radial position
    
    % Laminar flow profile (parabolic)
    V_avg = case_data.Velocity;
    V_laminar = 2 * V_avg * (1 - (r/R).^2);
    
    % Turbulent flow profile (1/7th power law)
    V_turbulent = V_avg * (7/6) * (1 - abs(r/R)).^(1/7);
    
    % Plot comparison
    subplot(2, 3, case_idx);
    plot(V_laminar * 100, r * 1000, 'b-', 'LineWidth', 2, 'DisplayName', 'Laminar');
    hold on;
    plot(V_turbulent * 100, r * 1000, 'r--', 'LineWidth', 2, 'DisplayName', 'Turbulent');
    xlabel('Velocity (cm/s)');
    ylabel('Radial Position (mm)');
    title(sprintf('%s\nD=%.1fmm, V=%.2fm/s, Re=%.0f', ...
        case_data.Fluid{1}, case_data.Diameter*1000, case_data.Velocity, case_data.Reynolds));
    legend('Location', 'best');
    grid on;
    
    % Plot normalized profiles
    subplot(2, 3, case_idx + 3);
    plot(V_laminar / max(V_laminar), r/R, 'b-', 'LineWidth', 2, 'DisplayName', 'Laminar');
    hold on;
    plot(V_turbulent / max(V_turbulent), r/R, 'r--', 'LineWidth', 2, 'DisplayName', 'Turbulent');
    xlabel('Normalized Velocity (V/V_{max})');
    ylabel('Normalized Position (r/R)');
    title('Normalized Velocity Profiles');
    legend('Location', 'best');
    grid on;
end

%% Part 3: Reynolds Number Effects

fprintf('\n=== PART 3: Reynolds Number Effects ===\n\n');

% Generate friction factor plot
Re_range = logspace(1, 8, 1000); % Reynolds number range

% Laminar flow friction factor
f_laminar = 64 ./ Re_range;

% Turbulent flow friction factor (Blasius)
f_turbulent = 0.316 ./ Re_range.^0.25;

% Colebrook-White for rough pipes (approximation)
f_colebrook = 0.11 * (68./Re_range + 0.0002).^0.25;

figure('Position', [100, 100, 1200, 500]);

% Friction factor plot
subplot(1, 2, 1);
loglog(Re_range, f_laminar, 'b-', 'LineWidth', 2, 'DisplayName', 'Laminar: f=64/Re');
hold on;
loglog(Re_range, f_turbulent, 'r-', 'LineWidth', 2, 'DisplayName', 'Turbulent (Blasius): f=0.316/Re^{0.25}');
loglog(Re_range, f_colebrook, 'g--', 'LineWidth', 2, 'DisplayName', 'Turbulent (Colebrook)');

% Mark transition region
transition_start = 2000;
transition_end = 4000;
x_fill = [transition_start, transition_end, transition_end, transition_start];
y_fill = [1e-3, 1e-3, 1e-1, 1e-1];
fill(x_fill, y_fill, 'yellow', 'FaceAlpha', 0.3, 'EdgeColor', 'none', ...
    'DisplayName', 'Transition Region');

xlabel('Reynolds Number (Re)');
ylabel('Darcy Friction Factor (f)');
title('Friction Factor vs Reynolds Number');
legend('Location', 'best');
grid on;
xlim([1e1, 1e8]);
ylim([1e-3, 1e-1]);

% Velocity profiles for different Re
subplot(1, 2, 2);
R = 0.05; % Pipe radius
r = linspace(-R, R, 100);
Re_values = [500, 2000, 5000, 10000, 50000];
colors = lines(length(Re_values));

for i = 1:length(Re_values)
    Re_val = Re_values(i);
    
    if Re_val < 2000
        % Laminar profile
        V_profile = 2 * (1 - (r/R).^2);
        style = '-';
    else
        % Turbulent profile - power law with exponent depending on Re
        if Re_val < 1e5
            n = 7; % 1/7 power law
        else
            n = 8 + 0.1*log10(Re_val/1e5); % Increasing exponent with Re
        end
        V_profile = ((n+1)/(2*n)) * (2)^(1/n) * (1 - abs(r/R)).^(1/n);
        style = '--';
    end
    
    plot(V_profile, r*1000, style, 'Color', colors(i,:), 'LineWidth', 2, ...
        'DisplayName', sprintf('Re = %d', Re_val));
    hold on;
end

xlabel('Normalized Velocity');
ylabel('Radial Position (mm)');
title('Velocity Profiles vs Reynolds Number');
legend('Location', 'best');
grid on;

%% Part 4: Interactive Calculator

fprintf('\n=== PART 4: Interactive Calculator ===\n\n');

% Create interactive GUI
callInteracrtiveCalculator();


%% Part 5: Flow Regime Concepts

fprintf('\n=== PART 5: Flow Regime Concepts ===\n\n');

% Create streamline plots for different flow regimes
createFlowVisualizations();

%% Bonus: CFD Simulation (Simplified 2D)

fprintf('\n=== BONUS: CFD Simulation ===\n\n');

% Perform simplified 2D CFD simulation
runCFDSimulation();

fprintf('\n=== ANALYSIS COMPLETE ===\n');

%% Supporting Functions

function callInteracrtiveCalculator()
    try
        system('node server.js &');  % run in background
        pause(2);  % wait for server to start
        web('http://localhost:3000/calculator.html', '-browser');  % adjust port if needed
    catch ME
        disp("Error launching server:");
        disp(ME.message);
    end
end

function createFlowVisualizations()
    % Create visualizations for different flow regimes
    
    figure('Position', [100, 100, 1200, 400]);
    sgtitle('Flow Regime Visualization', 'FontSize', 14, 'FontWeight', 'bold');
    
    % Laminar flow
    subplot(1, 3, 1);
    [X, Y] = meshgrid(0:0.1:10, -1:0.1:1);
    U = ones(size(X));
    V = zeros(size(X));
    streamslice(X, Y, U, V, 2);
    title('Laminar Flow');
    xlabel('Axial Direction');
    ylabel('Radial Direction');
    axis equal;
    
    % Transitional flow
    subplot(1, 3, 2);
    U_trans = ones(size(X)) + 0.3 * sin(2*pi*X/2) .* exp(-abs(Y));
    V_trans = 0.1 * cos(2*pi*X/2) .* Y;
    streamslice(X, Y, U_trans, V_trans, 2);
    title('Transitional Flow');
    xlabel('Axial Direction');
    ylabel('Radial Direction');
    axis equal;
    
    % Turbulent flow
    subplot(1, 3, 3);
    U_turb = ones(size(X));
    V_turb = zeros(size(X));
    
    % Add random fluctuations for turbulence
    for i = 1:size(X, 1)
        for j = 1:size(X, 2)
            fluctuation = 0.5 * randn * exp(-abs(Y(i,j)));
            U_turb(i,j) = U_turb(i,j) + fluctuation;
            if j > 1
                V_turb(i,j) = 0.2 * randn * exp(-abs(Y(i,j)));
            end
        end
    end
    
    streamslice(X, Y, U_turb, V_turb, 2);
    title('Turbulent Flow');
    xlabel('Axial Direction');
    ylabel('Radial Direction');
    axis equal;
end

function runCFDSimulation()
    % Simplified 2D CFD simulation for pipe flow
    
    fprintf('Running simplified 2D CFD simulation...\n');
    
    % Simulation parameters
    L = 10;        % Pipe length
    H = 2;         % Pipe height
    nx = 100;      % Grid points in x
    ny = 50;       % Grid points in y
    nu = 0.01;     % Kinematic viscosity
    U_inlet = 1;   % Inlet velocity
    
    % Initialize grid
    dx = L/(nx-1);
    dy = H/(ny-1);
    x = linspace(0, L, nx);
    y = linspace(-H/2, H/2, ny);
    [X, Y] = meshgrid(x, y);
    
    % Initialize velocity field
    U = zeros(ny, nx);
    V = zeros(ny, nx);
    
    % Boundary conditions
    U(:,1) = U_inlet * (1 - (2*Y(:,1)/H).^2); % Parabolic inlet profile
    U(:,end) = U(:,end-1); % Outlet
    U(1,:) = 0;   % No-slip at top wall
    U(end,:) = 0; % No-slip at bottom wall
    
    % Time parameters
    dt = 0.001;
    total_time = 5;
    n_steps = round(total_time/dt);
    
    % Simulation loop
    fprintf('Progress: ');
    for step = 1:n_steps
        if mod(step, n_steps/10) == 0
            fprintf('%d%% ', round(100*step/n_steps));
        end
        
        % Calculate intermediate velocity
        U_temp = U;
        
        % Convection and diffusion (simplified)
        for i = 2:ny-1
            for j = 2:nx-1
                % Simplified Navier-Stokes (viscous term only)
                d2u_dx2 = (U(i,j+1) - 2*U(i,j) + U(i,j-1))/dx^2;
                d2u_dy2 = (U(i+1,j) - 2*U(i,j) + U(i-1,j))/dy^2;
                
                U(i,j) = U_temp(i,j) + dt * nu * (d2u_dx2 + d2u_dy2);
            end
        end
        
        % Apply boundary conditions
        U(:,1) = U_inlet * (1 - (2*Y(:,1)/H).^2);
        U(:,end) = U(:,end-1);
        U(1,:) = 0;
        U(end,:) = 0;
    end
    fprintf('\n');
    
    % Plot results
    figure('Position', [100, 100, 800, 600]);
    
    subplot(2, 2, 1);
    contourf(X, Y, U, 20, 'LineColor', 'none');
    colorbar;
    title('Velocity Magnitude');
    xlabel('Axial Position');
    ylabel('Radial Position');
    
    subplot(2, 2, 2);
    streamslice(X, Y, U, V, 2);
    title('Streamlines');
    xlabel('Axial Position');
    ylabel('Radial Position');
    
    subplot(2, 2, 3);
    plot(y, U(:,round(nx/2)), 'b-', 'LineWidth', 2);
    hold on;
    plot(y, U(:,end), 'r--', 'LineWidth', 2);
    title('Velocity Profiles');
    xlabel('Radial Position');
    ylabel('Velocity');
    legend('Mid-pipe', 'Outlet');
    grid on;
    
    subplot(2, 2, 4);
    centerline_velocity = U(round(ny/2), :);
    plot(x, centerline_velocity, 'k-', 'LineWidth', 2);
    title('Centerline Velocity Development');
    xlabel('Axial Position');
    ylabel('Centerline Velocity');
    grid on;
    
    fprintf('CFD simulation completed successfully.\n');
end
