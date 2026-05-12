clc; 
clear variables; % Replaced clear all with clear variables for better performance
clf;
close all
% ==========================================
% USER OPTIONS
% ==========================================

% ANIMATION_MODE: 
% true  = Show real-time animation with 3D tube visualization during simulation
% false = Run simulation, then generate static subplot figures from new data
ANIMATION_MODE = true; 

% LOAD_OBSTACLES_FROM_MAT:
% true  = Load  obstacle definitions with 60 uncertain unsafe from 'mobile_robot_obstacles.mat'
% false = Generate random multiple oscillating obstacles dynamically
LOAD_OBSTACLES_FROM_MAT = true;

% PLOT_MODE: 
% 1 = Simple circular obstacles (no probabilistic regions)
% 2 = Probability contour obstacles (shows probabilistic regions \varepsilon^j,1.08\varepsilon^j sublevel sets)
PLOT_MODE = 1; 

% Snapshot times for static visualization (only used if ANIMATION_MODE = false)
snap_times = [3 6 13];
num_plots = length(snap_times);

% ==========================================
% INITIALIZATION
% ==========================================
cen = [0, 0]; 
robot = cen; % Start position
goal = [9.5, 9.5];  % Goal position

%% STT Parameters
k_1 = 0.5;   % goal reaching term constant
k_2 = 0.4;   % constant for probabilistic avoidance
k_3 = 2*k_2; % constant for probabilistic avoidance
rad_max = 0.8; % maximum tube radius
p2 = 0.9999; % probability within which obstacle avoidance should start
rad_min = 0.1; % minimum tube radius

%% Simulation Parameters
delt = 0.02; % Step size for movement
tf = 16;
tc = 0.9*tf;
t = 0:delt:tf;
n_o = 60; % Number of multiple random oscillating obstacles

N_map_colors = 100; % Number of steps in the colormap
light_red = [1 0.5 0.5]; % Lightest red/pink
dark_red = [1 0 0];   % Darkest red/maroon
obstacle_red_colormap = [ ...
    linspace(light_red(1), dark_red(1), N_map_colors)', ...
    linspace(light_red(2), dark_red(2), N_map_colors)', ...
    linspace(light_red(3), dark_red(3), N_map_colors)'  ...
];

% ==========================================
% OBSTACLE SETUP
% ==========================================
if LOAD_OBSTACLES_FROM_MAT
    filename = 'mobile_robot_obstacles.mat';
    if isfile(filename)
        data = load(filename);
        if isfield(data, 'obstacles_store')
            obstacles = data.obstacles_store;
            fprintf('60 uncertain obstacles successfully loaded from %s\n', filename);
            n_o = length(obstacles); % Ensure n_o matches loaded data
        else
            error('Variable obstacles_store not found in %s', filename);
        end
    else
        warning('File %s not found. Falling back to random obstacle generation.', filename);
        LOAD_OBSTACLES_FROM_MAT = false;
    end
end

if ~LOAD_OBSTACLES_FROM_MAT
    fprintf('Generating random oscillating obstacles...\n');
    rng('shuffle'); % Shuffle random seed for different results each run
    obstacles(n_o) = struct('mu', [], 'Sigma', [], 'p1', [], 'R', [],'v',[],'p_hist',[], ...
                            'Qmasked',[],'X1',[],'X2',[], 'obs_type', [], 'mu_center', [], 'amplitude', [], 'omega', []);
                        
    for i = 1:n_o
        valid_position = false;
        while ~valid_position
            % Random center position between x:[-3, 10] and y:[-3, 10]
            rand_mu = [rand()*(11+3) - 3, rand()*(11+3) - 3];
            % Ensure obstacles don't spawn too close to start or goal
            if norm(rand_mu - cen) > 2.5 && norm(rand_mu - goal) > 2.0
                valid_position = true;
            end
        end

        obstacles(i).obs_type = 'oscillating';
        obstacles(i).mu_center = rand_mu;

        % Random amplitude between [-2, 2] for both x and y
        amp_x = (rand() - 0.5) * 4; 
        amp_y = (rand() - 0.5) * 4;
        obstacles(i).amplitude = [amp_x, amp_y]; 

        % Random frequency between 0.4 and 1.2
        obstacles(i).omega = 0.4 + rand() * 0.8; 

        % Initial State
        obstacles(i).mu = obstacles(i).mu_center;
        obstacles(i).Sigma = 0.1 + rand() * (1.0-0.1); % Random Sigma [0.1, 1.0]
        obstacles(i).p1 = 0.7 + rand() * (0.99-0.7);  % Random probability [0.7, 0.99]
        obstacles(i).R = 0.2 + rand() *(0.5-0.2);     % Random radius [0.2, 0.5]
        obstacles(i).v = [0, 0]; 
        obstacles(i).p_hist = [];
    end
end

n = n_o; % Use n for loops

% ==========================================
% PRE-ALLOCATE TRAJECTORY ARRAYS
% ==========================================
num_steps = length(t);
robot_trajectory = zeros(num_steps, 3);
tube_trajectory = zeros(num_steps, 3);
radius_array = zeros(num_steps, 1);

% Initial states
robot_trajectory(1,:) = [robot(1), robot(2), t(1)];
tube_trajectory(1,:) = [cen(1), cen(2), t(1)];
radius_array(1) = rad_max;

% ==========================================
% ANIMATION PLOT SETUP (ONLY IF ANIMATION_MODE = TRUE)
% ==========================================
if ANIMATION_MODE
    disp("Starting real-time simulation and animation...")
    figure(1);
    set(gcf,"Position", [100 100 1200 600]);
    
    % --- 2D Plot Setup ---
    subplot(1,2,1)
    hold on; grid on; box on;
    set(gcf,'Color','k'); set(gca,'Color','k');
    axis([-4 12 -4 12]); axis square;
    ax = gca; ax.FontSize = 18;
    ax.XColor = 'w'; ax.YColor = 'w'; ax.GridColor = 'w'; ax.MinorGridColor = 'w'; ax.LineWidth = 1.5;
    xlabel('$x$ (m)','Interpreter','latex','Color','w');
    ylabel('$y$ (m)','Interpreter','latex','Color','w');
    
    % Draw Target
    theta_cir = linspace(0, 2*pi, 100);
    fill(goal(1) + 1.1*cos(theta_cir), goal(2) + 1.1*sin(theta_cir), 'g', 'FaceAlpha', 0.4, 'EdgeColor', 'g');
    
    % Init Obstacles
    obstacle_plots = cell(1, n);
    for o = 1:n
        if PLOT_MODE == 2
            %=========Pre-Computation of Probabilstoc Region==========
            % We will precompute the sub-level sets of \hatq which give us
            % the probabilistic regions of avoidance, this only works when uncertainitiny doesn't vary with time, in case of varying uncertainity go with Plot_Mode=1
            %uncomment the below line to generate proabbilstic region if running for different set of unsafe set other than present at working.mat file
            % [obstacles(o).Qmasked, obstacles(o).X1, obstacles(o).X2] = compute_q_region(obstacles(o).mu', obstacles(o).Sigma, obstacles(o).R, 1-obstacles(o).p1);
            [~, obstacle_plots{o}] = contourf( ...
                obstacles(o).X1, obstacles(o).X2, obstacles(o).Qmasked, ...
                [0  obstacles(o).p1 ], "ShowText", false, "FaceAlpha", 0.5); 
            colormap(gca, obstacle_red_colormap);
        else
            x_obs = obstacles(o).mu(1) + obstacles(o).R * cos(theta_cir);
            y_obs = obstacles(o).mu(2) + obstacles(o).R * sin(theta_cir);
            obstacle_plots{o} = fill(x_obs, y_obs, light_red, 'FaceAlpha', 0.5, 'EdgeColor', dark_red, 'LineWidth', 1.5);
        end
    end
    
    % Init Moving entities
    robot_plot = plot(robot(1), robot(2), 'co', 'MarkerFaceColor', 'c', 'MarkerSize', 3);
    tube_plot = plot(cen(1), cen(2), 'bo', 'MarkerSize', 3.8);
    trajectory_line = plot(robot(1), robot(2), 'b-', 'LineWidth', 1.5);
    
    % --- 3D Plot Setup ---
    subplot(1,2,2)
    hold on; grid on; box on;
    set(gca,'Color','k'); axis([-5 12 -5 12 0 tc]); axis square; view(3);
    ax2 = gca; ax2.FontSize = 18;
    ax2.XColor = 'w'; ax2.YColor = 'w'; ax2.ZColor = 'w'; ax2.GridColor = 'w'; ax2.MinorGridColor = 'w'; ax2.LineWidth = 1.5;
    xlabel('$x$ (m)','Interpreter','latex','Color','w');
    ylabel('$y$ (m)','Interpreter','latex','Color','w');
    zlabel('$t$ (s)','Interpreter','latex','Color','w');
    
    % 3D Target
    [target_x, target_y, target_z] = cylinder(1, 50);
    target_z = target_z * 1 + (tc - 1);
    surf(target_x+goal(1), target_y+goal(2), target_z, 'FaceColor', 'g', 'FaceAlpha', 0.4, 'EdgeColor', 'none');
    
    robot_path_handle = plot3(robot_trajectory(1,1), robot_trajectory(1,2), robot_trajectory(1,3), 'c-', 'LineWidth', 2);
    tube_center_handle = plot3(tube_trajectory(1,1), tube_trajectory(1,2), tube_trajectory(1,3), 'b-', 'LineWidth', 1);
    
    [X_cyl, Y_cyl, Z_cyl] = cylinder(rad_max, 50);
    safety_tube_handle = surf(X_cyl + cen(1), Y_cyl + cen(2), Z_cyl * delt + t(1), ...
        'FaceColor',[0.3010 0.7450 0.9330], 'FaceAlpha', 0.4, 'EdgeColor', 'none');
else
    disp("Running simulation for static visualization...")
end

% ==========================================
% MAIN SIMULATION LOOP (ALWAYS RUNS)
% ==========================================
stop_idx = num_steps;

for iter = 1:num_steps
    current_time = t(iter);
    
    d_cen_1 = -k_1 * sign(cen - goal).*abs(cen - goal).^(1/3); %% target driven terms
    d_cen_2 = [0, 0];
    rad = rad_max; 
    sum_dist = 10^4; 
    
    %% probabilistic avoidance term calculations
    for o = 1:n
        mu = obstacles(o).mu;
        s = obstacles(o).Sigma;
        varepsilon = obstacles(o).p1;
        lambda = (norm(mu-cen)/s)^2;
        Rs = obstacles(o).R + rad_min;
        p_i = ncx2cdf((Rs/s)^2, 2, lambda);
        q = 1 - p_i;
        delta = cen - mu;
        
        if q > varepsilon && q < p2
            m_j = (1/q - 1/p2) * (1 / (q-varepsilon)) * ( delta')';
            v_j = [0 1; -1 0] * m_j';
            d_cen_2 = d_cen_2 + (k_2*m_j + k_3*(v_j)');
        end
        quant = s*sqrt((ncx2inv(1-varepsilon, 2, lambda))) - obstacles(o).R;
        sum_dist = min(sum_dist, quant);
    end
    
    rad = min(rad_max, sum_dist);
    radius_array(iter) = rad;
    
    % Update probabilities and kinematics for next step
    for o = 1:n
        % Calculate probability update
        mu = obstacles(o).mu;
        s = obstacles(o).Sigma;
        varepsilon = obstacles(o).p1;
        lambda = (norm(mu-cen)/s)^2;
        Rs = obstacles(o).R + rad_min;
        p_update = 1 - ncx2cdf(((Rs-rad_min)/s)^2, 2, lambda);
        if p_update < varepsilon
            disp("failure")
        end
        obstacles(o).p_hist = [obstacles(o).p_hist p_update];
        
        % Update obstacle position dynamically
        obstacles(o).mu = obstacles(o).mu_center + obstacles(o).amplitude * sin(obstacles(o).omega * current_time);
        obstacles(o).v = obstacles(o).amplitude * obstacles(o).omega * cos(obstacles(o).omega * current_time);
    end
    
    d_cen_total = d_cen_1 + d_cen_2;
    
    % Update positions
    if norm(d_cen_total) > 0
        cen = cen + delt * (d_cen_total);
    end
    
    % Replace control() with your actual robot control logic. 
    % Assuming standard pursuit logic for demonstration if control() is a custom function.
    u_robot = control(cen, rad, robot); 
    robot = robot + delt * u_robot;
    
    % Save state
    robot_trajectory(iter,:) = [robot(1), robot(2), current_time];
    tube_trajectory(iter,:) = [cen(1), cen(2), current_time];
    
    % --- ANIMATION UPDATE ---
    if ANIMATION_MODE
        subplot(1,2,1);
        
        for o=1:n
            if PLOT_MODE == 2
                x1 = obstacles(o).X1 + obstacles(o).v(1,1)*delt;
                y2 = obstacles(o).X2 + obstacles(o).v(1,2)*delt;
                obstacles(o).X1 = x1;
                obstacles(o).X2 = y2;
                
                if ishandle(obstacle_plots{o}), delete(obstacle_plots{o}); end
                Q1 = obstacles(o).Qmasked;
                varepsilon = obstacles(o).p1;
                Q1(Q1 > min(1.08*varepsilon,0.995)) = NaN;    
                [~, obstacle_plots{o}] = contourf(obstacles(o).X1, obstacles(o).X2, Q1, ...
                    [0 varepsilon min(1.08*varepsilon,0.995)], 'FaceAlpha',0.5, 'LineStyle','none'); 
                colormap(gca, obstacle_red_colormap);
            else
                x_obs = obstacles(o).mu(1) + obstacles(o).R * cos(theta_cir);
                y_obs = obstacles(o).mu(2) + obstacles(o).R * sin(theta_cir);
                if ishandle(obstacle_plots{o}) && isvalid(obstacle_plots{o})
                    set(obstacle_plots{o}, 'XData', x_obs, 'YData', y_obs);
                else
                    obstacle_plots{o} = fill(x_obs, y_obs, light_red, 'FaceAlpha', 0.5, 'EdgeColor', dark_red, 'LineWidth', 1.5);
                end
            end
        end
        
        set(robot_plot, 'XData', robot(1), 'YData', robot(2));
        set(tube_plot, 'XData', cen(1), 'YData', cen(2));
        set(trajectory_line, 'XData', robot_trajectory(1:iter,1), 'YData', robot_trajectory(1:iter,2));
        
        delete(findobj('Type', 'line', 'Color', [0.3010 0.7450 0.9330]));
        plot(cen(1) + (rad - 0.1)*cos(theta_cir), cen(2) + (rad - 0.1)*sin(theta_cir), ...
            'color', [0.3010 0.7450 0.9330], 'LineStyle', '--', 'LineWidth', 1.5);
            
        subplot(1,2,2);
        set(robot_path_handle, 'XData', robot_trajectory(1:iter,1), 'YData', robot_trajectory(1:iter,2), 'ZData', robot_trajectory(1:iter,3));
        set(tube_center_handle, 'XData', tube_trajectory(1:iter,1), 'YData', tube_trajectory(1:iter,2), 'ZData', tube_trajectory(1:iter,3));
        
        Z_scaled = Z_cyl * delt + current_time;
        set(safety_tube_handle, 'XData', X_cyl*rad + cen(1), 'YData', Y_cyl*rad + cen(2), 'ZData', Z_scaled);
        
        if iter > 1
            surf(X_cyl*rad + cen(1), Y_cyl*rad + cen(2), Z_scaled, 'FaceColor', [0.3010 0.7450 0.9330], 'FaceAlpha', 0.4, 'EdgeColor', 'none');
        end
        
        [az, el] = view; view(az + 120/tf*delt, el);
        drawnow;
    end
    
    % Goal Check
    if norm(cen - goal) < 0.001
        disp(['Goal reached at t = ', num2str(current_time)]);
        stop_idx = iter;
        break;
    end
end

% Trim arrays to actual ending iteration length
t = t(1:stop_idx);
robot_trajectory = robot_trajectory(1:stop_idx, :);
tube_trajectory = tube_trajectory(1:stop_idx, :);
radius_array = radius_array(1:stop_idx);

% ==========================================
% STATIC SUBPLOT VISUALIZATION (IF NOT ANIMATED)
% ==========================================
if ~ANIMATION_MODE
    disp("Generating static subplot visualization from newly calculated simulation...")
    r = 1.1;
    figure(2);
    set(gcf,"Color",'w', "Position", [100 100 1400 400]);
    theta = linspace(0, 2*pi, 100);
    
    for i = 1:num_plots
        [~, idx] = min(abs(t - snap_times(i)));
        snap_tube = tube_trajectory(idx,1:2);
        snap_rad  = radius_array(idx);
        
        subplot(1, num_plots+1, i)
        hold on; grid on;
        axis([-4 12 -4 12]); 
        title(sprintf('Time = %.2f sec', snap_times(i)));
        xlabel('$x$(m)','Interpreter','latex'); 
        ylabel('$y$(m)','Interpreter','latex');
        
        % Target
        fill(goal(1) + r*cos(theta), goal(2) + r*sin(theta), 'g','FaceAlpha',0.43,'EdgeColor','none');
        
        % Tube safety circle
        fill(snap_tube(1) + (snap_rad - 0.1)*cos(theta), snap_tube(2) + (snap_rad - 0.1)*sin(theta), ...
            'c', 'FaceAlpha',0.54,'EdgeColor','c')
        
        % Trajectory until this time
        plot(robot_trajectory(1:idx,1), robot_trajectory(1:idx,2), 'k', 'LineWidth',1.5);
        plot(tube_trajectory(1:idx,1), tube_trajectory(1:idx,2), 'b', 'LineWidth',1.5); 
        plot(robot_trajectory(idx,1), robot_trajectory(idx,2), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 4);
        
        % Obstacles at snap_time
        for o = 1:n_o
            obs = obstacles(o);
            mu = obs.mu_center + obs.amplitude .* sin(obs.omega * snap_times(i));
            
            if PLOT_MODE == 2
                Q = obs.Qmasked;
                varepsilon = obs.p1;
                Q1 = Q;
                Q1(Q > min(1.08*varepsilon,0.995)) = NaN;
                
                % Recalculate grid shift for snapshot time (assuming initial grid centered at mu_center)
                X1_snap = obs.X1 + (mu(1) - obs.mu_center(1));
                X2_snap = obs.X2 + (mu(2) - obs.mu_center(2));
                
                contourf(X1_snap, X2_snap, Q1, [0 varepsilon min(1.08*varepsilon,0.995)], 'FaceAlpha',0.5, 'LineStyle','none'); 
                colormap(gca, obstacle_red_colormap);
            else
                fill(mu(1) + obs.R * cos(theta), mu(2) + obs.R * sin(theta), ...
                    light_red, 'FaceAlpha', 0.5, 'EdgeColor', dark_red, 'LineWidth', 1.5);
            end
        end
    end
    
    % 3D Spatiotemporal Tube Subplot
    subplot(1,num_plots+1,num_plots+1)
    hold on; grid on; view(3);
    title('3D Spatiotemporal Tube','FontSize',14);
    xlabel('$x$(m)','Interpreter','latex'); 
    ylabel('$y$(m)','Interpreter','latex'); 
    zlabel('$t$(s)','Interpreter','latex');
    
    [Xc, Yc, Zc] = cylinder(1, 30); 
    for i = 1:3:length(t) % Tube skip
        R_tube = radius_array(i);
        surf(Xc * R_tube + tube_trajectory(i,1), Yc * R_tube + tube_trajectory(i,2), Zc * delt + tube_trajectory(i,3), ...
            'FaceColor', [0 0.4 1], 'FaceAlpha', 0.5, 'EdgeColor','none');
    end
    plot3(robot_trajectory(:,1), robot_trajectory(:,2), t, 'Color', 'k', 'LineWidth', 1.7);
    plot3(tube_trajectory(:,1), tube_trajectory(:,2), t, 'Color', 'b', 'LineWidth', 1.7);
end

% ==========================================
% PROBABILITY HISTORY PLOT (ALWAYS RUNS)
% ==========================================
figure(3);
hold on; grid on;
colors = lines(n_o); 
for i = 1:n_o
    % Match history length to actual iteration count
    hist_len = min(length(t), length(obstacles(i).p_hist));
    plot(t(1:hist_len), obstacles(i).p_hist(1:hist_len)', 'LineWidth', 1.5, 'Color', colors(i,:), ...
        'DisplayName', ['Obs ', num2str(i), ' p(t)']);
    yline(obstacles(i).p1, '--', 'Color', colors(i,:), 'LineWidth', 1, 'HandleVisibility', 'off'); 
end
title('Probability Evolution for All Obstacles');
xlabel('Time (s)');
ylabel('Probability p(t)');
legend('show', 'Location', 'best');
hold off;