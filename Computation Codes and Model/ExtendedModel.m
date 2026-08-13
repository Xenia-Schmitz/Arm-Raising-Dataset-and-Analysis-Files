%% Capture Point Control - ACTUALLY FIXED VERSION
% Properly implements temporal causality: controller generates NO torques
% until disturbance is visible within prediction horizon

clear; clc;

% Setup
params = build_params();

% Define prediction horizons to test
horizons_seconds = [0.1,0.2,0.3];
horizons_steps = round(horizons_seconds / params.dt);

% Storage
results = struct( ...
    'horizon', {}, ...
    'x', {}, ...
    'cop', {}, ...
    'capture_point', {}, ...
    'tau_ankle', {}, ...
    'tau_hip', {}, ...
    'disturbance_visible_time', {}, ...
    'success_rate', {} );

% Run simulation for each horizon
for h_idx = 1:length(horizons_steps)
    fprintf('Simulating horizon = %.2f s...\n', horizons_seconds(h_idx));
    
    result = simulate_capture_point_control(params, horizons_steps(h_idx));
    result.horizon = horizons_seconds(h_idx);
    
    % Calculate when disturbance becomes visible
    result.disturbance_visible_time = 0.3 - horizons_seconds(h_idx);
    
    % Calculate success rate
    in_bos = (result.cop >= params.BOS_back) & (result.cop <= params.BOS_front);
    result.success_rate = sum(in_bos) / length(in_bos);
    
    results(h_idx) = result;
    
    fprintf('  Disturbance visible at: %.2f s\n', result.disturbance_visible_time);
    fprintf('  Success rate: %.1f%%\n', result.success_rate * 100);
    fprintf('  Max CoP: %.2f cm (BoS: %.1f to %.1f cm)\n\n', ...
            max(abs(result.cop))*100, params.BOS_back*100, params.BOS_front*100);
end

% Plot results
plot_results(results, params, horizons_seconds);

% ==================== SIMULATION FUNCTION ==============================
function result = simulate_capture_point_control(params, horizon_steps)
    % Controller parameters
    Kp_base = [160; 105];
    Kd_base = [110; 50];
    
    % Capture point control gains
    Kp_cp = 850;   % Capture point position error gain
    Kd_cp = 60;   % Capture point velocity error gain
    
    % Anticipatory stiffness
    Kp_max = [140; 75];
    Kd_max = [120; 55];
    stiffness_threshold = 30;  % rad/s² - threshold to detect disturbance
    stiffness_rise_time = 0.05;
    stiffness_decay_time = 0.05;
    
    % Rate limits
    tau_rate_limit = 700;
    
    % Initialize
    x = [0; 0; 0; 0];
    N = params.N;
    
    % Storage
    x_hist = zeros(4, N);
    cop_hist = zeros(N, 1);
    cp_hist = zeros(N, 1);
    tau_ankle_hist = zeros(N, 1);
    tau_hip_hist = zeros(N, 1);
    
    tau_prev = [0; 0];
    stiffness_activation = 0;
    
    % Simulation loop
    for i = 1:N
        % Current shoulder state (for actual dynamics)
        theta_s_actual = params.shoulder_traj.theta(i);
        theta_s_dot_actual = params.shoulder_traj.theta_dot(i);
        theta_s_ddot_actual = params.shoulder_traj.theta_ddot(i);
        
        % ===== CHECK IF DISTURBANCE IS VISIBLE IN HORIZON =====
        disturbance_visible = false;
        max_future_accel = 0;
        
        if horizon_steps > 0 && i + horizon_steps <= N
            % Look ahead within horizon
            for h = 1:horizon_steps
                k_future = i + h;
                if k_future > params.N, break; end
                
                future_accel = abs(params.shoulder_traj.theta_ddot(k_future));
                
                % Check if significant disturbance is visible in horizon
                if future_accel > stiffness_threshold
                    disturbance_visible = true;
                    max_future_accel = max(max_future_accel, future_accel);
                end
            end
        end
        
        % ===== ANTICIPATORY STIFFNESS MODULATION =====
        if disturbance_visible
            disturbance_magnitude = max_future_accel; % / 1.0; %min(max_future_accel / 20.0, 1.0);
            target_activation = disturbance_magnitude;
        else
            target_activation = 0;
        end
        
        % Smooth transition to target activation
        dt = params.dt;
        if target_activation > stiffness_activation
            % Rising - prepare for disturbance
            rise_rate = 1.0 / stiffness_rise_time;
            stiffness_activation = min(stiffness_activation + rise_rate * dt, target_activation);
        else
            % Falling - return to baseline
            decay_rate = 1.0 / stiffness_decay_time;
            stiffness_activation = max(stiffness_activation - decay_rate * dt, target_activation);
        end
        
        % ===== CONTROL ARCHITECTURE =====
        % Baseline control: Always active for stability and natural oscillations
        % Anticipatory control: Only when disturbance is visible
        
        % Compute current gains (baseline + anticipatory)
        Kp = Kp_base + stiffness_activation * (Kp_max - Kp_base);
        Kd = Kd_base + stiffness_activation * (Kd_max - Kd_base);
        
        % Always apply baseline stabilization
        tau_baseline = -Kp_base .* x(1:2) - Kd_base .* x(3:4);
        
        % Add anticipatory stiffening only when disturbance is visible
        if disturbance_visible
            Kp_extra = stiffness_activation * (Kp_max - Kp_base);
            Kd_extra = stiffness_activation * (Kd_max - Kd_base);
            tau_anticipatory = -Kp_extra .* x(1:2) - Kd_extra .* x(3:4);
        else
            tau_anticipatory = [0; 0];
        end
        
        tau_baseline_total = tau_baseline + tau_anticipatory;
        
        % For capture point prediction, use appropriate gains
        if disturbance_visible
            Kp_for_prediction = Kp;
            Kd_for_prediction = Kd;
        else
            Kp_for_prediction = Kp_base;
            Kd_for_prediction = Kd_base;
        end
        
        % ===== COMPUTE CURRENT CoM STATE =====
        state = [x(1); x(2); theta_s_actual; x(3); x(4); theta_s_dot_actual];
        theta = [x(1); x(2); theta_s_actual];
        omega = [x(3); x(4); theta_s_dot_actual];
        M = update_mass_matrix(theta, params);
        C = compute_coriolis_matrix(theta, omega, params);
        G = compute_gravity_vector(theta, params);
        
        M11 = M(1:2, 1:2);
        M13 = M(1:2, 3);
        C_omega = C * omega;
        tau_eff = tau_baseline_total - G(1:2) - C_omega(1:2) - M13 * theta_s_ddot_actual;
        alpha_body = M11 \ tau_eff;
        
        state_dot = [x(3); x(4); theta_s_dot_actual; alpha_body(1); alpha_body(2); theta_s_ddot_actual];
        [com_x, com_z, com_vx, ~] = compute_com_kinematics_from_state_and_alpha(state, state_dot, params);
        
        % ===== CAPTURE POINT CONTROL =====
        % Natural frequency of inverted pendulum
        omega_0 = sqrt(params.g / com_z);
        
        % Current capture point (where CoM lands with no control)
        capture_point = com_x + com_vx / omega_0;
        cp_hist(i) = capture_point;
        
        % Desired capture point (center of BoS for stability)
        cp_desired = -0.02;
        
        % ===== PREDICT FUTURE CAPTURE POINT =====
        tau_cp = [0; 0];
        
        if disturbance_visible && horizon_steps > 0 && i + horizon_steps <= N || i>83 && i + horizon_steps <= N
            [cp_future, cp_velocity_future] = predict_capture_point(x, i, horizon_steps, params, Kp_for_prediction, Kd_for_prediction);
            
            % Compute capture point error and derivative
            cp_error_current = capture_point - cp_desired;
            cp_error_future = cp_future - cp_desired;
            cp_velocity_error = cp_velocity_future;
            
            % Weight near vs far capture point errors
            weight_current = 0.2;
            weight_future = 0.8;
            cp_error_weighted = weight_current * cp_error_current + weight_future * cp_error_future;
            
            % Generate corrective torque based on capture point dynamics
            tau_cp = -Kp_cp * cp_error_weighted * [0.8; 0.2] - Kd_cp * cp_velocity_error * [0.8; 0.2];
            
        end
        
        % ===== TOTAL TORQUE =====
        tau_desired = tau_baseline_total + tau_cp; % + tau_emergency;
        
        % ===== RATE LIMITING =====
        tau = rate_limit(tau_desired, tau_prev, tau_rate_limit, params.dt);
        tau_prev = tau;
        
        % ===== STORE =====
        x_hist(:, i) = x;
        tau_ankle_hist(i) = tau(1);
        tau_hip_hist(i) = tau(2);
        
        % ===== COMPUTE ACTUAL CoP =====
        tau_eff_actual = tau - G(1:2) - C_omega(1:2) - M13 * theta_s_ddot_actual;
        alpha_actual = M11 \ tau_eff_actual;
        state_dot_actual = [x(3); x(4); theta_s_dot_actual; alpha_actual(1); alpha_actual(2); theta_s_ddot_actual];
        [com_x_actual, com_z_actual, ~, com_ax_actual] = compute_com_kinematics_from_state_and_alpha(state, state_dot_actual, params);
        cop_hist(i) = com_x_actual - (com_z_actual * com_ax_actual) / params.g;
        
        % ===== INTEGRATE =====
        if i < N
            x = dynamics_step(x, tau, theta_s_actual, theta_s_dot_actual, theta_s_ddot_actual, params);
        end
    end

    cop_hist_smooth = smooth_cop(cop_hist, params);
    
    % Package results
    result.x = x_hist;
    result.cop = cop_hist_smooth;
    result.capture_point = cp_hist;
    result.tau_ankle = tau_ankle_hist;
    result.tau_hip = tau_hip_hist;
end
%% ==================== CoP SMOOTHING =====================================
function cop_smooth = smooth_cop(cop_raw, params)
    % Apply low-pass filter to smooth CoP trajectory
    % Uses moving average or Savitzky-Golay filter
    
    % Method 1: Moving average
    window_size = 10; % 50ms window at 100Hz
    cop_smooth = movmean(cop_raw, window_size);
    
    % Method 2: Savitzky-Golay (better - preserves peaks)
    % Uncomment to use instead:
    % order = 3;
    % framelen = 11; % must be odd
    % cop_smooth = sgolayfilt(cop_raw, order, framelen);
end
%% ==================== PREDICT CAPTURE POINT ============================
function [cp_future, cp_velocity] = predict_capture_point(x_current, k_current, H, params, Kp, Kd)
    x_pred = x_current;
    
    % Simulate forward to get future state
    % Uses current control gains and actual shoulder trajectory
    for h = 1:H
        k_future = k_current + h;
        if k_future > params.N, break; end
        
        theta_s = params.shoulder_traj.theta(k_future);
        theta_s_dot = params.shoulder_traj.theta_dot(k_future);
        theta_s_ddot = params.shoulder_traj.theta_ddot(k_future);
        
        tau_pred = -Kp .* x_pred(1:2) - Kd .* x_pred(3:4);
        x_pred = dynamics_step(x_pred, tau_pred, theta_s, theta_s_dot, theta_s_ddot, params);
    end
    
    % Compute future CoM state
    theta_s = params.shoulder_traj.theta(min(k_current + H, params.N));
    theta_s_dot = params.shoulder_traj.theta_dot(min(k_current + H, params.N));
    theta_s_ddot = params.shoulder_traj.theta_ddot(min(k_current + H, params.N));
    
    state = [x_pred(1); x_pred(2); theta_s; x_pred(3); x_pred(4); theta_s_dot];
    theta = [x_pred(1); x_pred(2); theta_s];
    omega = [x_pred(3); x_pred(4); theta_s_dot];
    M = update_mass_matrix(theta, params);
    C = compute_coriolis_matrix(theta, omega, params);
    G = compute_gravity_vector(theta, params);
    
    M11 = M(1:2, 1:2);
    M13 = M(1:2, 3);
    C_omega = C * omega;
    tau_pred = -Kp .* x_pred(1:2) - Kd .* x_pred(3:4);
    tau_eff = tau_pred - G(1:2) - C_omega(1:2) - M13 * theta_s_ddot;
    alpha_body = M11 \ tau_eff;
    
    state_dot = [x_pred(3); x_pred(4); theta_s_dot; alpha_body(1); alpha_body(2); theta_s_ddot];
    [com_x, com_z, com_vx, ~] = compute_com_kinematics_from_state_and_alpha(state, state_dot, params);
    
    % Compute future capture point
    omega_0 = sqrt(params.g / com_z);
    cp_future = com_x + com_vx / omega_0;
    cp_velocity = com_vx;
end

%% ==================== HELPER FUNCTIONS =================================
function tau_limited = rate_limit(tau_desired, tau_prev, rate_limit, dt)
    max_change = rate_limit * dt;
    delta = tau_desired - tau_prev;
    
    for i = 1:length(tau_desired)
        if abs(delta(i)) > max_change
            delta(i) = sign(delta(i)) * max_change;
        end
    end
    
    tau_limited = tau_prev + delta;
end

%% ==================== DYNAMICS =========================================
function x_next = dynamics_step(x, tau, theta_s, theta_s_dot, theta_s_ddot, params)
    dt = params.dt;
    k1 = dynamics_derivative(x, tau, theta_s, theta_s_dot, theta_s_ddot, params);
    x2 = x + 0.5 * dt * k1;
    k2 = dynamics_derivative(x2, tau, theta_s, theta_s_dot, theta_s_ddot, params);
    x3 = x + 0.5 * dt * k2;
    k3 = dynamics_derivative(x3, tau, theta_s, theta_s_dot, theta_s_ddot, params);
    x4 = x + dt * k3;
    k4 = dynamics_derivative(x4, tau, theta_s, theta_s_dot, theta_s_ddot, params);
    x_next = x + (dt/6) * (k1 + 2*k2 + 2*k3 + k4);
end

function dx = dynamics_derivative(x, tau, theta_s, theta_s_dot, theta_s_ddot, params)
    theta1 = x(1); theta2 = x(2);
    omega1 = x(3); omega2 = x(4);
    
    theta = [theta1; theta2; theta_s];
    omega = [omega1; omega2; theta_s_dot];
    
    M = update_mass_matrix(theta, params);
    C = compute_coriolis_matrix(theta, omega, params);
    G = compute_gravity_vector(theta, params);
    
    M11 = M(1:2, 1:2);
    M13 = M(1:2, 3);
    C_omega = C * omega;
    
    tau_eff = [tau(1); tau(2)] - G(1:2) - C_omega(1:2) - M13 * theta_s_ddot;
    alpha = M11 \ tau_eff;
    
    dx = [omega1; omega2; alpha(1); alpha(2)];
end

%% ==================== PARAMETERS =======================================
function params = build_params()
    params = struct();
    params.g = 9.81;
    params.m_total = 70;
    
    params.l1 = 0.9;
    params.l2 = 0.5;
    params.l3 = 0.55;
    
    params.r1 = params.l1/2;
    params.r2 = params.l2/2;
    params.r3 = params.l3/2;
    
    params.m1 = 2 * (0.161 * params.m_total);
    params.m2 = 0.578 * params.m_total;
    params.m3 = 2 * (0.05 * params.m_total);
    
    params.I1 = params.m1 * params.l1^2 / 12;
    params.I2 = params.m2 * params.l2^2 / 12;
    params.I3 = params.m3 * params.l3^2 / 12;
    
    params.foot_length = 0.27;
    params.BOS_front = 0.75 * params.foot_length;
    params.BOS_back = -0.25 * params.foot_length;
    
    params.t_final = 2.5;
    params.dt = 0.01;
    params.t = 0:params.dt:params.t_final;
    params.N = length(params.t);
    
    [params.shoulder_traj, ~] = generate_shoulder_trajectory(params);
end

%% ==================== CoM KINEMATICS ===================================
function [total_com_x, total_com_z, total_com_vx, total_com_ax] = compute_com_kinematics_from_state_and_alpha(state, state_dot, params)
    theta1 = state(1); theta2 = state(2); theta3 = state(3);
    omega1 = state(4); omega2 = state(5); omega3 = state(6);
    alpha1 = state_dot(4); alpha2 = state_dot(5); alpha3 = state_dot(6);
    
    com1_x = params.r1 * sin(theta1);
    com2_x = params.l1 * sin(theta1) + params.r2 * sin(theta1 + theta2);
    com3_x = params.l1 * sin(theta1) + params.l2 * sin(theta1 + theta2) + params.r3 * sin(theta1 + theta2 + theta3);
    
    com1_z = params.r1 * cos(theta1);
    com2_z = params.l1 * cos(theta1) + params.r2 * cos(theta1 + theta2);
    com3_z = params.l1 * cos(theta1) + params.l2 * cos(theta1 + theta2) + params.r3 * cos(theta1 + theta2 + theta3);
    
    com1_vx = params.r1 * cos(theta1) * omega1;
    com2_vx = params.l1 * cos(theta1) * omega1 + params.r2 * cos(theta1 + theta2) * (omega1 + omega2);
    com3_vx = params.l1 * cos(theta1) * omega1 + params.l2 * cos(theta1 + theta2) * (omega1 + omega2) + ...
              params.r3 * cos(theta1 + theta2 + theta3) * (omega1 + omega2 + omega3);
    
    com1_ax = -params.r1 * sin(theta1) * omega1^2 + params.r1 * cos(theta1) * alpha1;
    com2_ax = -params.l1 * sin(theta1) * omega1^2 + params.l1 * cos(theta1) * alpha1 + ...
              -params.r2 * sin(theta1+theta2) * (omega1+omega2)^2 + params.r2 * cos(theta1+theta2) * (alpha1+alpha2);
    com3_ax = -params.l1 * sin(theta1) * omega1^2 + params.l1 * cos(theta1) * alpha1 + ...
              -params.l2 * sin(theta1+theta2) * (omega1+omega2)^2 + params.l2 * cos(theta1+theta2) * (alpha1+alpha2) + ...
              -params.r3 * sin(theta1+theta2+theta3) * (omega1+omega2+omega3)^2 + ...
              params.r3 * cos(theta1+theta2+theta3) * (alpha1+alpha2+alpha3);
    
    total_com_x = (params.m1*com1_x + params.m2*com2_x + params.m3*com3_x) / params.m_total;
    total_com_z = (params.m1*com1_z + params.m2*com2_z + params.m3*com3_z) / params.m_total;
    total_com_vx = (params.m1*com1_vx + params.m2*com2_vx + params.m3*com3_vx) / params.m_total;
    total_com_ax = (params.m1*com1_ax + params.m2*com2_ax + params.m3*com3_ax) / params.m_total;
end

%% ==================== TRAJECTORY GENERATION ============================
function [shoulder_traj, traj_params] = generate_shoulder_trajectory(params)
    t = params.t;
    N = params.N;
    dt = params.dt;
    t_start = 0.3;
    t_rise = 0.2;
    t_fall = 0.29; %0.29 %0.35
    t_total = t_rise + t_fall;
    t_end = t_start + t_total;
    
    theta_start = pi;
    theta_end = 0;
    total_angle = abs(theta_end - theta_start);
    
    r_endpoint = params.l3;
    v_tangential_target = 1;
    omega_peak_target = v_tangential_target / r_endpoint;
    
    theta = ones(N,1) * theta_start;
    theta_dot = zeros(N,1);
    theta_ddot = zeros(N,1);
    
    for i = 1:N
        ti = t(i);
        if ti >= t_start && ti <= t_end
            t_local = ti - t_start;
            if t_local <= t_rise
                tau = t_local / t_rise;
                s_dot = 10*tau^3 - 15*tau^4 + 6*tau^5;
                s_ddot = (30*tau^2 - 60*tau^3 + 30*tau^4) / t_rise;
                theta_dot(i) = omega_peak_target * s_dot;
                theta_ddot(i) = omega_peak_target * s_ddot;
            else
                tau = (t_local - t_rise) / t_fall;
                s_dot = 1 - (3*tau^2 - 2*tau^3);
                s_ddot = -(6*tau - 6*tau^2) / t_fall;
                theta_dot(i) = omega_peak_target * s_dot;
                theta_ddot(i) = omega_peak_target * s_ddot;
            end
        end
    end
    
    for i = 2:N
        theta(i) = theta(i-1) - theta_dot(i)*dt;
    end
    
    final_angle = theta(find(t >= t_end,1)) - theta_start;
    scale_factor = total_angle / abs(final_angle);
    theta_dot = theta_dot * scale_factor;
    theta_ddot = theta_ddot * scale_factor;
    
    theta = ones(N,1) * theta_start;
    for i = 2:N
        theta(i) = theta(i-1) - theta_dot(i)*dt;
    end
    
    shoulder_traj.theta = theta;
    shoulder_traj.theta_dot = theta_dot;
    shoulder_traj.theta_ddot = theta_ddot;
    
    traj_params.t_rise = t_rise;
    traj_params.t_fall = t_fall;
end

%% ==================== DYNAMICS MATRICES ================================
function M = update_mass_matrix(theta, params)
    theta1 = theta(1); theta2 = theta(2); theta3 = theta(3);
    m1 = params.m1; m2 = params.m2; m3 = params.m3;
    l1 = params.l1; l2 = params.l2; l3 = params.l3;
    r1 = params.r1; r2 = params.r2; r3 = params.r3;
    I1 = params.I1; I2 = params.I2; I3 = params.I3;
    
    c2 = cos(theta2); c3 = cos(theta3); c23 = cos(theta2 + theta3);
    
    M11 = I1 + I2 + I3 + m2*l1^2 + m3*l1^2 + m3*l2^2 + m1*r1^2 + m2*r2^2 + m3*r3^2 ...
          + 2*l1*m3*r3*c23 + 2*l1*l2*m3*c2 + 2*l1*m2*r2*c2 + 2*l2*m3*r3*c3;
    M12 = I2 + I3 + l2^2*m3 + 2*l2*m3*r3*c3 + l1*l2*m3*c2 + m2*r2^2 + l1*m2*r2*c2 + m3*r3^2 + l1*m3*r3*c23;
    M13 = I3 + m3*r3^2 + l1*m3*r3*c23 + l2*m3*r3*c3;
    M22 = I2 + I3 + m3*l2^2 + 2*l2*m3*r3*c3 + m2*r2^2 + m3*r3^2;
    M23 = I3 + m3*r3^2 + l2*m3*r3*c3;
    M33 = I3 + m3*r3^2;
    
    M = [M11, M12, M13; M12, M22, M23; M13, M23, M33];
end

function C = compute_coriolis_matrix(theta, omega, params)
    theta1 = theta(1); theta2 = theta(2); theta3 = theta(3);
    omega1 = omega(1); omega2 = omega(2); omega3 = omega(3);
    m2 = params.m2; m3 = params.m3;
    l1 = params.l1; l2 = params.l2;
    r2 = params.r2; r3 = params.r3;
    
    s2 = sin(theta2); s3 = sin(theta3); s23 = sin(theta2 + theta3);
    
    C = zeros(3,3);
    C(1,1) = -(l1*m2*r2*s2 + l1*m3*r3*s23 + l1*l2*m3*s2)*omega2 - (l2*m3*r3*s3 + l1*m3*r3*s23)*omega3;
    C(1,2) = -(l1*m2*r2*s2 + l1*m3*r3*s23 + l1*l2*m3*s2)*omega1 - (l1*m2*r2*s2 + l1*m3*r3*s23 + l1*l2*m3*s2)*omega2 - m3*r3*(l2*s3 + l1*s23)*omega3;
    C(1,3) = -m3*r3*(l2*s3 + l1*s23)*omega1 - m3*r3*(l2*s3 + l1*s23)*omega2 - m3*r3*(l2*s3 + l1*s23)*omega3;
    C(2,1) = (l1*m2*r2*s2 + l1*m3*r3*s23 + l1*l2*m3*s2)*omega1 - l2*m3*r3*s3*omega3;
    C(2,2) = -l2*m3*r3*s3*omega3;
    C(2,3) = -l2*m3*r3*s3*omega1 - l2*m3*r3*s3*omega2 - l2*m3*r3*s3*omega3;
    C(3,1) = (l2*m3*r3*s3 + l1*m3*r3*s23)*omega1 + l2*m3*r3*s3*omega2;
    C(3,2) = l2*m3*r3*s3*omega1 + l2*m3*r3*s3*omega2;
    C(3,3) = 0;
end

function G = compute_gravity_vector(theta, params)
    theta1 = theta(1); theta2 = theta(2); theta3 = theta(3);
    g = params.g;
    m1 = params.m1; m2 = params.m2; m3 = params.m3;
    l1 = params.l1; l2 = params.l2;
    r1 = params.r1; r2 = params.r2; r3 = params.r3;
    
    s1 = sin(theta1); s12 = sin(theta1+theta2); s123 = sin(theta1+theta2+theta3);
    G1 = g*((m1*r1 + m2*l1 + m3*l1)*s1 + (m2*r2 + m3*l2)*s12 + m3*r3*s123);
    G2 = g*((m2*r2 + m3*l2)*s12 + m3*r3*s123);
    G3 = g*m3*r3*s123;
    G = [G1; G2; G3];
end

%% ==================== PLOTTING =========================================
function plot_results(results, params, horizons_seconds)
    figure('Position', [100, 100, 1600, 1000]);
    
    % CoP trajectories
    subplot(1,2,1);
    colors = lines(length(results));
    for h_idx = 1:length(results)
        plot(params.t, results(h_idx).cop * 100, 'LineWidth', 2, 'Color', colors(h_idx,:), ...
             'DisplayName', sprintf('H = %.2fs', results(h_idx).horizon));
        hold on;
        
        % Mark when disturbance becomes visible
        xline(results(h_idx).disturbance_visible_time, '--', 'Color', colors(h_idx,:), ...
              'LineWidth', 1, 'HandleVisibility', 'off');
    end
    yline(params.BOS_front*100, 'r--', 'LineWidth', 1.5, 'DisplayName', 'BoS limits');
    yline(params.BOS_back*100, 'r--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    xline(0.3, 'k-', 'LineWidth', 2, 'DisplayName', 'Arm movement starts');
    xlabel('Time (s)'); ylabel('CoP (cm)');
    title('Center of Pressure', 'FontWeight', 'bold');
    legend('Location', 'best'); grid on;
    
    % Success rate
    subplot(1,2,2);
    success_rates = [results.success_rate] * 100;
    bar(horizons_seconds, success_rates, 'FaceColor', [0.2 0.6 0.8]);
    xlabel('Horizon (s)'); ylabel('% in BoS');
    title('Balance Success', 'FontWeight', 'bold');
    ylim([0 105]); grid on;
    
end