% Clear workspace and set random seed for reproducibility
clear all;
clc;
rng(123); % Set seed for reproducibility

% ### Step 1: Define Simulation Parameters (Table I from Paper)
N = 16;              % Number of BS antennas
K_n = 5;             % Number of known interferers in the same cell
K_u_vec = [20, 30, 40]; % Number of unknown interferers in neighboring cells
r_vec = [100, 200];  % Distances between desired UE and BS (m)
alpha = 3.76;        % Pathloss exponent (not used in new model)
P_tx = 0.1;          % Transmit power per UE (100 mW = 0.1 W)
noise_power = 10^(-8.5) * 10^(-3); % -85 dBm in W (-94 dBm + 9 dB noise figure)
tau_c = 200;         % Coherence block length
tau_p = 10;          % Pilot sequence length
tau_u = tau_c - tau_p; % Data transmission length
num_realizations = 1000; % Number of Monte Carlo realizations
fc = 2e9;            % Carrier frequency (2 GHz)

% ### Step 2: Spatial Layout
theta_known = 2 * pi * rand(1, K_n); % Random angles for known UEs
r_known = [60, 100, 140, 180, 220];  % Radii for known UEs (m)
x_known = r_known .* cos(theta_known);
y_known = r_known .* sin(theta_known);
d_known = sqrt(x_known.^2 + y_known.^2);

% ### Step 3: Precompute Pathloss for Known Interferers
PL_known_dB = 20 + 35 * log10(d_known); % Adjusted pathloss model
PL_known_linear = 10.^(-PL_known_dB / 10); % Large-scale fading coefficients (beta_i)

% ### Step 4: Monte Carlo Simulation
SINR_RZF = zeros(length(r_vec), length(K_u_vec), num_realizations);
SINR_MR = zeros(length(r_vec), length(K_u_vec), num_realizations); % For MR combining
SINR_b_RZF = zeros(length(r_vec), length(K_u_vec), num_realizations); % For baseline
SE_RZF = zeros(length(r_vec), length(K_u_vec), num_realizations); % Store SE for plotting

for r_idx = 1:length(r_vec)
    r = r_vec(r_idx);
    % Pathloss for desired UE
    PL_desired_dB = 20 + 35 * log10(r);
    PL_desired_linear = 10^(-PL_desired_dB / 10);
    
    for k_idx = 1:length(K_u_vec)
        K_u = K_u_vec(k_idx);
        % Generate unknown interferers
        r_unknown = 250 + (500 - 250) * rand(1, K_u);
        theta_unknown = 2 * pi * rand(1, K_u);
        x_unknown = r_unknown .* cos(theta_unknown);
        y_unknown = r_unknown .* sin(theta_unknown);
        d_unknown = sqrt(x_unknown.^2 + y_unknown.^2);
        PL_unknown_dB = 20 + 35 * log10(d_unknown);
        PL_unknown_linear = 10.^(-PL_unknown_dB / 10);
        
        for realization = 1:num_realizations
            % Generate true channels
            h_desired = sqrt(PL_desired_linear) * (randn(N, 1) + 1j*randn(N, 1))/sqrt(2);
            h_known = sqrt(PL_known_linear) .* (randn(N, K_n) + 1j*randn(N, K_n))/sqrt(2);
            h_unknown = sqrt(PL_unknown_linear) .* (randn(N, K_u) + 1j*randn(N, K_u))/sqrt(2);
            
            % Channel Estimation (Eq. 2)
            pilot_interferers = rand(1, K_u) < (tau_p / (K_n + K_u + 1)); % Random pilot assignment
            h_pilot_interf = h_unknown(:, pilot_interferers);
            y_pilot = sqrt(tau_p * P_tx) * h_desired + ...
                      sqrt(tau_p * P_tx) * sum(h_pilot_interf, 2) + ...
                      sqrt(noise_power) * (randn(N, 1) + 1j*randn(N, 1))/sqrt(2);
            R_k = PL_desired_linear * eye(N); % Covariance matrix = beta_k * I_N
            Psi = tau_p * P_tx * R_k + noise_power * eye(N);
            h_hat_desired = sqrt(tau_p * P_tx) * R_k * (Psi \ y_pilot);
            
            % Channel estimates for known interferers (assuming orthogonal pilots)
            h_hat_known = zeros(N, K_n);
            for i = 1:K_n
                y_pilot_i = sqrt(tau_p * P_tx) * h_known(:,i) + ...
                            sqrt(noise_power) * (randn(N, 1) + 1j*randn(N, 1))/sqrt(2);
                R_i = PL_known_linear(i) * eye(N);
                Psi_i = tau_p * P_tx * R_i + noise_power * eye(N);
                h_hat_known(:,i) = sqrt(tau_p * P_tx) * R_i * (Psi_i \ y_pilot_i);
            end
            
            % RZF Combining (Eq. 10)
            h_hat_all = [h_hat_desired, h_hat_known];
            R = P_tx * (h_hat_all * h_hat_all') + noise_power * eye(N);
            v_RZF = (R \ (P_tx * h_hat_desired));
            v_RZF = v_RZF / norm(v_RZF); % Normalize
            
            % MR Combining
            v_MR = h_hat_desired / norm(h_hat_desired); % Normalize
            
            % SINR Calculation for RZF (Eq. 6)
            DS_k = sqrt(P_tx) * (v_RZF' * h_desired); % Desired signal
            IUI_u = P_tx * sum(abs(v_RZF' * h_unknown).^2); % Unknown interference
            h_all_known = [h_desired, h_known]; % Set D_n
            IUSI_u = P_tx * sum(abs(v_RZF' * h_all_known).^2) - P_tx * abs(v_RZF' * h_desired)^2;
            sigma_tilde = noise_power * norm(v_RZF)^2;
            SINR_RZF(r_idx, k_idx, realization) = abs(DS_k)^2 / (IUI_u + IUSI_u + sigma_tilde);
            
            % SINR Calculation for MR
            DS_k_MR = sqrt(P_tx) * (v_MR' * h_desired);
            IUI_u_MR = P_tx * sum(abs(v_MR' * h_unknown).^2);
            IUSI_u_MR = P_tx * sum(abs(v_MR' * h_all_known).^2) - P_tx * abs(v_MR' * h_desired)^2;
            sigma_tilde_MR = noise_power * norm(v_MR)^2;
            SINR_MR(r_idx, k_idx, realization) = abs(DS_k_MR)^2 / (IUI_u_MR + IUSI_u_MR + sigma_tilde_MR);
            
            % Baseline SINR (Eq. 19)
            SINR_b_RZF(r_idx, k_idx, realization) = abs(DS_k)^2 / (IUSI_u + sigma_tilde);
            
            % Spectral Efficiency (Eq. 5)
            SE_RZF(r_idx, k_idx, realization) = (tau_u / tau_c) * log2(1 + SINR_RZF(r_idx, k_idx, realization));
        end
    end
end

% Compute ε-outage SE 
epsilon_values = 0.05:0.05:0.3;
num_eps = length(epsilon_values);
SE_proposed = zeros(length(r_vec), length(K_u_vec), num_eps);
m_values = [3.10, 3.42]; % Margin values from paper
SE_baseline = zeros(length(r_vec), length(K_u_vec), length(m_values));
epsilon_baseline = zeros(length(r_vec), length(K_u_vec), length(m_values));

% Proposed method
for r_idx = 1:length(r_vec)
    for k_idx = 1:length(K_u_vec)
        sorted_SE = sort(squeeze(SE_RZF(r_idx, k_idx, :)));
        for e_idx = 1:num_eps
            T_idx = ceil(epsilon_values(e_idx) * num_realizations);
            SE_proposed(r_idx, k_idx, e_idx) = sorted_SE(max(1, T_idx));
        end
    end
end

% Baseline method
for m_idx = 1:length(m_values)
    m = m_values(m_idx);
    for r_idx = 1:length(r_vec)
        for k_idx = 1:length(K_u_vec)
            avg_SINR_b = mean(squeeze(SINR_b_RZF(r_idx, k_idx, :)));
            SE_b = (tau_u/tau_c) * log2(1 + avg_SINR_b/m);
            SE_baseline(r_idx, k_idx, m_idx) = SE_b;
            outage = mean(squeeze(SINR_RZF(r_idx, k_idx, :)) < avg_SINR_b/m);
            epsilon_baseline(r_idx, k_idx, m_idx) = outage;
        end
    end
end

% Plot
figure;
hold on;
markers = {'o', 's', 'd'};
colors_r = lines(length(r_vec));
colors_m = {'r', 'b'}; % Different colors for m=3.10 and m=3.42

% Plot proposed curves
for r_idx = 1:length(r_vec)
    for k_idx = 1:length(K_u_vec)
        plot(epsilon_values, squeeze(SE_proposed(r_idx, k_idx, :)), ...
            'Color', colors_r(r_idx,:), 'LineStyle', '-', ...
            'Marker', markers{k_idx}, 'LineWidth', 2);
    end
end

% Plot baseline markers
for m_idx = 1:length(m_values)
    for r_idx = 1:length(r_vec)
        for k_idx = 1:length(K_u_vec)
            eps = epsilon_baseline(r_idx, k_idx, m_idx);
            if eps <= 0.3 % Exclude points with epsilon > 0.3
                scatter(eps, SE_baseline(r_idx, k_idx, m_idx), ...
                    100, colors_m{m_idx}, markers{k_idx}, 'filled', 'LineWidth', 2);
            end
        end
    end
end

xlabel('Outage probability (\epsilon)');
ylabel('Spectral Efficiency (bps/Hz)');
title('\epsilon-outage SE with RZF Combining');
%legendEntries = cell(1, 6 + 12); % 6 proposed + 12 baseline (2 margins * 6)

legendEntries = {
    'Proposed K_u=20, r=100m'
    'Proposed K_u=30, r=100m'
    'Proposed K_u=40, r=100m'
    'Proposed K_u=20, r=200m'
    'Proposed K_u=40, r=200m'
    'Proposed K_u=30, r=200m'
    'Baseline m=3.10, K_u=20, r=100m'
    'Baseline m=3.10, K_u=30, r=100m'
    'Baseline m=3.10, K_u=40, r=100m'
    'Baseline m=3.10, K_u=20, r=200m'
    'Baseline m=3.10, K_u=30, r=200m'
    'Baseline m=3.10, K_u=40, r=200m'
    'Baseline m=3.42, K_u=20, r=100m'
    'Baseline m=3.42, K_u=30, r=100m'
    'Baseline m=3.42, K_u=40, r=100m'
    'Baseline m=3.42, K_u=20, r=200m'
    'Baseline m=3.42, K_u=30, r=200m'
    'Baseline m=3.42, K_u=40, r=200m'
};
legend(legendEntries, 'Location', 'bestoutside', 'NumColumns', 2);
grid on;
hold off;