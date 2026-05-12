%% STT + Dynamic Obstacles in 3D
clc;
clear all;
clf;
% close all

% ==========================================================
% VISUALIZATION MODE
% ==========================================================
% true  -> Animation
% false -> Static subplot visualization
ANIMATION_MODE = false;

% ==========================================
% STATIC PLOT OPTIONS
% ==========================================
snap_times = [ 5 10 34];
num_plots = length(snap_times);

%% --- Parameters ---
cen = [1, 1, 1];
robot = cen;
velo = [0 0 0];
goal = [8, 8, 8];
%% tim evarying bounds parameters
rho_L = -2;
rho_U = 2;
rho_d = rho_U - rho_L;

%% STT Pramaters
k2 = 0.03;
k3 = k2;
k1 = 0.31;
rad_max = 0.9;
p2_d = 0.99999;
rad_min = 0.3;

%% Simulation Parameters
delt = 0.009;
tf = 41;
tc = 0.9 * tf;
t = 0:delt:tf;
N = 4;
time = [];
% ==========================================================
% OBSTACLES Information
% ==========================================================
obstacles(N) = struct('mu', [], 'Sigma', [], 'p1', [], ...
    'R', [], 'v', [], 'p_hist', [], 'c', []);

obstacles(1).mu = [10, 10, 10];
obstacles(1).Sigma = 1;
obstacles(1).p1 = 0.99;
obstacles(1).R = 0.5;
obstacles(1).v = [-0.11, -0.11, -0.11];
obstacles(1).p_hist = [];
obstacles(1).c = [1.0, 0.0, 0.0];

obstacles(2).mu = [3, 3, 3];
obstacles(2).Sigma = 0.02;
obstacles(2).p1 = 0.7;
obstacles(2).R = 1;
obstacles(2).v = [-0.2, -0.2, -0.2];
obstacles(2).p_hist = [];
obstacles(2).c = [0.93, 0.44, 0.44];

obstacles(3).mu = [8, 8, 8];
obstacles(3).Sigma = 0.02;
obstacles(3).p1 = 0.7;
obstacles(3).R = 0.5;
obstacles(3).v = [-0.2, -0.2, -0.2];
obstacles(3).p_hist = [];
obstacles(3).c = [0.93, 0.44, 0.44];

obstacles(4).mu = [5.2, 5.2, 0];
obstacles(4).Sigma = 0.06;
obstacles(4).p1 = 0.9;
obstacles(4).R = 1.5;
obstacles(4).v = [0.0, 0.0, 0.4];
obstacles(4).p_hist = [];
obstacles(4).c = [0.9, 0.25, 0.25];

% copy for plotting
obstacles_2 = obstacles;

% ==========================================================
% PLOT SETUP
% ==========================================================
if ANIMATION_MODE
    figure(1)
    hold on;
    grid on;
    box on;
    axis([0 10 0 10 0 10]);
    axis equal
    set(gca, ...
        'XLimMode', 'manual', ...
        'YLimMode', 'manual', ...
        'ZLimMode', 'manual');
    ax = gca;
    ax.FontSize = 14;
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    view(3);
    % Goal Plotting
    [sx, sy, sz] = sphere(20);
    surf(goal(1) + sx, ...
         goal(2) + sy, ...
         goal(3) + sz, ...
         'FaceAlpha', 0.4, ...
         'EdgeColor', 'none', ...
         'FaceColor', 'g');
    % Obstacles
    obstacle_plots = gobjects(size(obstacles,1),1);
    for o = 1:N
        [sx, sy, sz] = sphere(20);
        x = obstacles(o).mu(1) + obstacles(o).R * sx;
        y = obstacles(o).mu(2) + obstacles(o).R * sy;
        z = obstacles(o).mu(3) + obstacles(o).R * sz;
        obstacle_plots(o) = surf(x, y, z, ...
            'FaceAlpha', 0.55, ...
            'EdgeColor', 'none', ...
            'FaceColor', obstacles(o).c);
    end
    % Robot plot
    robot_plot = plot3(robot(1), robot(2), robot(3), ...
        'co', 'MarkerFaceColor', 'c');
    tube_plot = plot3(cen(1), cen(2), cen(3), 'bo');
    trajectory = plot3(robot(1), robot(2), robot(3), ...
        'b-', 'LineWidth', 1.5);
    trajectory_tube = plot3(cen(1), cen(2), cen(3), ...
        'c-', 'LineWidth', 1.5);
end
tube_trajectory = [];
robot_trajectory = [];
rad_array = [];
%% ==========================================================
%% STT LOOP
%% ==========================================================
for iter = 1:length(t)
    tic  
    % Update obstacle positions
    for o = 1:N
        obstacles(o).mu = obstacles(o).mu + obstacles(o).v * delt;
        if ANIMATION_MODE
            [sx, sy, sz] = sphere(20);
            x = obstacles(o).mu(1) + obstacles(o).R * sx;
            y = obstacles(o).mu(2) + obstacles(o).R * sy;
            z = obstacles(o).mu(3) + obstacles(o).R * sz;
            set(obstacle_plots(o), ...
                'XData', x, ...
                'YData', y, ...
                'ZData', z);
        end
    end
    % ======================================================
    % Goal reaching term
    % ======================================================
    d_cen1 = -k1 * sign(cen - goal) .* abs(cen - goal).^(1/3);
    % ======================================================
    % Obstacle avoidance term
    % ======================================================
    d_cen2 = [0, 0, 0];
    rad = rad_max;
    sum_dist = 1e4;
    for o = 1:N
        mu = obstacles(o).mu;
        s = obstacles(o).Sigma;
        p1 = obstacles(o).p1;
        lambda = (norm(mu - cen))^2 / s;
        Rs = obstacles(o).R + rad_min;
        p_i = ncx2cdf((Rs)^2 / s, 3, lambda);
        q = 1 - p_i;
        delta = cen - mu;
        if q > p1 && q < p2_d
            m_j = (1/q - 1/p2_d) * k2 * ...
                (1 / (q - p1)) * (delta')';
            Null = null(m_j);
            if ~isempty(Null)
                v_j = k3 * Null(:,1)';
            else
                v_j = [0 0 0];
            end
            d_cen2 = d_cen2 + m_j + v_j;
        end
        quant = sqrt(s) * ...
            sqrt((ncx2inv(1 - p1, 3, lambda))) ...
            - obstacles(o).R;
        sum_dist = min(sum_dist, quant);
    end
    sim_time = toc;
    rad = min(rad_max, sum_dist);
    for o = 1:N
        mu_p = obstacles(o).mu;
        s_p = obstacles(o).Sigma;
        lambda_p = (norm(mu_p - cen))^2 / s_p;
        Rs_p = obstacles(o).R + rad_min;
        p_update = 1 - ncx2cdf( ...
            ((Rs_p - rad_min + rad))^2 / s_p, ...
            3, ...
            lambda_p);
        obstacles(o).p_hist = ...
            [obstacles(o).p_hist p_update];
    end
    d_cen_total = d_cen1 + d_cen2;

    if norm(d_cen_total) > 0
        cen = cen + delt * d_cen_total;
    end
    rho_L = -3 * exp(-0.03 * t(iter)) * ones(1,3);
    rho_U = 3 * exp(-0.03 * t(iter)) * ones(1,3);
    rho_d = (rho_U - rho_L);
    rho_d_dot = ...
        -3 * 0.03 * 2 * exp(-0.03 * t(iter)) ...
        * ones(1,3);
    % ======================================================
    % Control update
    % ======================================================
    u = control_3d( ...
        cen, ...
        rad, ...
        robot, ...
        velo, ...
        rho_L, ...
        rho_U, ...
        rho_d, ...
        rho_d_dot);

    robot = robot + delt * velo;
    velo = velo + delt * u;
    % ======================================================
    % Store trajectories
    % ======================================================
    time = [time t(iter)];
    tube_trajectory = [tube_trajectory; cen];
    rad_array = [rad_array rad];
    robot_trajectory = [robot_trajectory; robot];
    % ======================================================
    % Animation updates
    % ======================================================
    if ANIMATION_MODE
        set(robot_plot, ...
            'XData', robot(1), ...
            'YData', robot(2), ...
            'ZData', robot(3));
        set(tube_plot, ...
            'XData', cen(1), ...
            'YData', cen(2), ...
            'ZData', cen(3));
        set(trajectory, ...
            'XData', [get(trajectory,'XData'), robot(1)], ...
            'YData', [get(trajectory,'YData'), robot(2)], ...
            'ZData', [get(trajectory,'ZData'), robot(3)]);
        set(trajectory_tube, ...
            'XData', [get(trajectory_tube,'XData'), cen(1)], ...
            'YData', [get(trajectory_tube,'YData'), cen(2)], ...
            'ZData', [get(trajectory_tube,'ZData'), cen(3)]);
        % Safety sphere
        delete(findobj(gca, 'Tag', 'safe_sphere'));
        [sx, sy, sz] = sphere(10);
        surf(cen(1) + (rad - 0.1) * sx, ...
             cen(2) + (rad - 0.1) * sy, ...
             cen(3) + (rad - 0.1) * sz, ...
             'FaceAlpha', 0.4, ...
             'EdgeColor', 'none', ...
             'FaceColor', [0.3010 0.7450 0.9330], ...
             'Tag', 'safe_sphere');
        drawnow;
        [az, el] = view;
        view(az - 120/tc*0.02, el);
        pause(0.0005);
    end
    % ======================================================
    % Goal reached
    % ======================================================
    if norm(cen - goal) < 0.2
        disp(['Goal reached at t = ', num2str(t(iter))]);
        break;
    end
end

%% ==========================================================
%% STATIC SUBPLOT VISUALIZATION
%% ==========================================================
if ~ANIMATION_MODE
    figure
    for i = 1:num_plots
        [~, idx] = min(abs(t - snap_times(i)));
        snap_tube = tube_trajectory(idx,1:3);
        snap_robot = robot_trajectory(idx,1:3);
        snap_rad = rad_array(idx);
        subplot(1, num_plots, i)
        hold on;
        grid on;
        box on;
        axis([0 10 0 10 0 10]);
        axis equal
        xlabel('$x$(m)', 'Interpreter', 'latex');
        ylabel('$y$(m)', 'Interpreter', 'latex');
        zlabel('$z$(m)', 'Interpreter', 'latex');
        % Goal
        [sx, sy, sz] = sphere(20);
        surf(goal(1) + sx, ...
             goal(2) + sy, ...
             goal(3) + sz, ...
             'FaceAlpha', 0.2, ...
             'EdgeColor', 'none', ...
             'FaceColor', 'g');
        % Robot trajectory
        plot3(robot_trajectory(1:idx,1), ...
              robot_trajectory(1:idx,2), ...
              robot_trajectory(1:idx,3), ...
              'k', 'LineWidth', 1.5);
        % Tube trajectory
        plot3(tube_trajectory(1:idx,1), ...
              tube_trajectory(1:idx,2), ...
              tube_trajectory(1:idx,3), ...
              'b', 'LineWidth', 1.5);
        % Robot point
        plot3(snap_robot(1), ...
              snap_robot(2), ...
              snap_robot(3), ...
              'ko', ...
              'MarkerFaceColor', 'k');
        % Tube sphere
        [sx, sy, sz] = sphere(20);
        surf(snap_tube(1) + snap_rad * sx, ...
             snap_tube(2) + snap_rad * sy, ...
             snap_tube(3) + snap_rad * sz, ...
             'FaceAlpha', 0.2, ...
             'EdgeColor', 'none', ...
             'FaceColor', [0.3010 0.7450 0.9330]);
        % Obstacles
        for o = 1:N
            obs_mu = ...
                obstacles_2(o).mu + ...
                obstacles_2(o).v * snap_times(i);
            [sx, sy, sz] = sphere(20);
            x = obs_mu(1) + obstacles(o).R * sx;
            y = obs_mu(2) + obstacles(o).R * sy;
            z = obs_mu(3) + obstacles(o).R * sz;
            surf(x, y, z, ...
                'FaceAlpha', 0.5, ...
                'EdgeColor', 'none', ...
                'FaceColor', obstacles(o).c);
        end
        view(3)
        hold off
    end
   
end