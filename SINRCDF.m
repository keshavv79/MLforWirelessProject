% Clear workspace and set random seed for reproducibility
clear all;
clc;
rng(123); % Set seed for reproducibility

% ### Step 1: Define Simulation Parameters
N = 16;              % Number of BS antennas
K_n = 5;             % Number of known interferers in the same cell
K_u_vec = [20, 30, 40]; % Number of unknown interferers in neighboring cells
r_vec = [100, 200];  % Distances between desired UE and BS (m)
alpha = 3.76;        % Pathloss exponent
P_tx = 0.1;          % Transmit power per UE (W)
noise_power = 10^(-9.4); % Noise power (-94 dBm in W)
sigma2 = noise_power;
tau_c = 200;         % Coherence block length
tau_p = 10;          % Pilot sequence length
tau_u = tau_c - tau_p; % Data transmission length
num_realizations = 3000; % Number of Monte Carlo realizations
epsilon_values = linspace(0.01, 0.3, 20); % Outage probability (x-axis)
m = 3.42; % Fixed fade margin for baseline


% ### Step 2: Spatial Layout
theta_known = 2 * pi * rand(1, K_n); % Random angles for known UEs
r_known = [60, 100, 140, 180, 220];  % Radii for known UEs (m)
x_known = r_known .* cos(theta_known);
y_known = r_known .* sin(theta_known);
d_known = sqrt(x_known.^2 + y_known.^2);

% ### Step 3: Precompute Pathloss for Known Interferers
fc = 2e9; % Carrier frequency (Hz)
shadow_std_dev = 4;  % in dB
F_kl = shadow_std_dev * randn(1, K_n);  % Shadow fading per interferer
PL_known_dB = -30.5 - 36.7 * log10(d_known) + F_kl/10;
PL_known_linear = 10.^(-PL_known_dB / 10);


% ### Step 4: Monte Carlo Simulation
SINR_RZF = zeros(length(r_vec), length(K_u_vec), num_realizations); %2x3x1000
SINR_MR = zeros(length(r_vec), length(K_u_vec), num_realizations); %2x3x1000
SINR_b_RZF = zeros(length(r_vec), length(K_u_vec), num_realizations); % For baseline
epsilonOutage = zeros(length(r_vec), length(K_u_vec), num_realizations);
% Shadowing parameter
for r_idx = 1:length(r_vec)
    r = r_vec(r_idx);
    
    % Pathloss for desired UE with shadow fading
    F_desired = shadow_std_dev * randn();  % Single UE
    PL_desired_dB = -30.5 - 36.7 * log10(r) + F_desired/10;
    PL_desired_linear = 10^(-PL_desired_dB / 10);

    for k_idx = 1:length(K_u_vec)
        K_u = K_u_vec(k_idx);
        SE_proposed = zeros(size(epsilon_values));
        SE_baseline = zeros(size(epsilon_values));
   
        for realization = 1:num_realizations
            % Unknown interferers: random placement
            r_unknown = 250 + (500 - 250) * rand(1, K_u);
            theta_unknown = 2 * pi * rand(1, K_u);
            x_unknown = r_unknown .* cos(theta_unknown);
            y_unknown = r_unknown .* sin(theta_unknown);
            d_unknown = sqrt(x_unknown.^2 + y_unknown.^2);

            % Pathloss with shadow fading for unknown interferers
            F_unk = shadow_std_dev * randn(1, K_u);
            PL_unknown_dB = -30.5 - 36.7 * log10(d_unknown) + F_unk/10;
            PL_unknown_linear = 10.^(-PL_unknown_dB / 10);

            % Generate channels
            %h_desired = sqrt(PL_desired_linear) * (randn(N, 1) + 1j * randn(N, 1)) / sqrt(2);
            %h_known = sqrt(PL_known_linear) .* (randn(N, K_n) + 1j * randn(N, K_n)) / sqrt(2);
            % Generate unknown interferer channels
            h_unknown = (randn(N, K_u) + 1j*randn(N, K_u)) / sqrt(2);
            h_unknown = h_unknown .* sqrt(reshape(PL_unknown_linear, [1, K_u])); 

            % === Step 1: True channels (Rayleigh fading with pathloss) ===
            h_desired_true = sqrt(PL_desired_linear) * (randn(N,1) + 1j*randn(N,1)) / sqrt(2);
            h_known_true = sqrt(PL_known_linear) .* (randn(N, K_n) + 1j*randn(N, K_n)) / sqrt(2);
            
            % === Step 2: MMSE channel estimation (i.i.d. Rayleigh) ===
            % Pilot observation from same-cell UEs only (orthogonal pilots assumed)
            %noise_desired = sqrt(sigma2 / tau_p) * (randn(N,1) + 1j*randn(N,1)) / sqrt(2);
            %noise_known = sqrt(sigma2 / tau_p) * (randn(N, K_n) + 1j*randn(N, K_n)) / sqrt(2);
            
            yp_desired = sqrt(tau_p * P_tx) * h_desired_true + noise_power;
            yp_known = sqrt(tau_p * P_tx) * h_known_true + noise_power;
            
            % MMSE estimator with R_k = I
            scaling = sqrt(tau_p * P_tx) / (tau_p * P_tx + sigma2);
            h_desired = scaling * yp_desired;       % Estimated desired channel
            h_known = scaling * yp_known;           % Estimated known channels

            
            % RZF combining
            R = P_tx*(h_known_true * h_known_true') +eye(N);
            v_RZF = P_tx*(R \ h_desired_true) ;
            v_RZF = v_RZF / norm(v_RZF);

            v_MR = h_desired/((norm(h_desired_true))^2);
            v_MR = v_MR / norm(v_MR);
            
            % SINR Calculation
            desired_power = abs(v_RZF' * h_desired_true)^2;
            interf_known =   sum(P_tx*abs(v_RZF' * h_known_true).^2);
            interf_unknown =  (P_tx^2.5)*(sum(abs(v_RZF' * h_unknown).^2));
            noise = noise_power * (norm(v_RZF)^2);
            SINR_RZF(r_idx, k_idx, realization) = desired_power / (interf_known + interf_unknown + noise);
            
            desired_power_MR = P_tx*abs(v_MR' * h_desired)^2;
            interf_known_MR = P_tx * sum(abs(v_MR' * h_known).^2);
            interf_unknown_MR = P_tx * sum(abs(v_MR' * h_unknown).^2);
            noise_MR = noise_power * norm(v_MR)^2;
            SINR_MR(r_idx, k_idx, realization) = desired_power_MR / (interf_known_MR  + noise_MR);

            if realization == 1 && r_idx == 1 && k_idx == 1
            fprintf("desired = %.3e, known = %.3e, unknown = %.3e, noise = %.3e\n", ...
                desired_power, interf_known, interf_unknown, noise);
            end
        end
        
    % Plot proposed scheme
   end
end

% ### Step 5: Plot CDFs (Figure 2)
disp('Min/Max SINR in dB for each setting:');
for r_idx = 1:length(r_vec)
    for k_idx = 1:length(K_u_vec)
        SINR_data = squeeze(SINR_RZF(r_idx, k_idx, :));
        SINR_dB = 10 * log10(SINR_data);
        fprintf('r = %d, K_u = %d → min = %.2f dB, max = %.2f dB\n', ...
            r_vec(r_idx), K_u_vec(k_idx), min(SINR_dB), max(SINR_dB));
    end
end

figure;
hold on;
colors = lines(length(K_u_vec));
linestyles = {'-', '--'};
legendEntries = {};

for r_idx = 1:length(r_vec)
    for k_idx = 1:length(K_u_vec)
        % === RZF ===
        SINR_data = squeeze(SINR_RZF(r_idx, k_idx, :));
        SINR_dB = 10 * log10(SINR_data); % Convert to dB
        [F, X] = ecdf(SINR_dB);          % ECDF

        plot(X, F, 'Color', colors(k_idx,:), 'LineStyle', linestyles{r_idx}, 'LineWidth', 2);
        legendEntries{end+1} = sprintf('RZF: K_u=%d, r=%dm', K_u_vec(k_idx), r_vec(r_idx));

        % === MR ===
        SINR_data_MR = squeeze(SINR_MR(r_idx, k_idx, :));
        SINR_dB_MR = 10 * log10(SINR_data_MR);
        [F_MR, X_MR] = ecdf(SINR_dB_MR);
        if r_idx == 1 && k_idx == 1
            plot(X_MR, F_MR, 'Color', colors(k_idx,:), 'LineStyle', linestyles{r_idx}, 'LineWidth', 1, 'Marker', 'none', 'LineStyle', ':');
            legendEntries{end+1} = sprintf('MR: K_u=%d, r=%dm', K_u_vec(k_idx), r_vec(r_idx));
        end
    end
end

legend(legendEntries, 'Location', 'southeast', 'NumColumns', 2);
xlabel('SINR (dB)');
ylabel('CDF');
title('CDF of SINR (dB) for RZF and MR Combining');
grid on;
set(gca, 'XScale', 'linear'); % dB is already logarithmic

% Overlay Analytical CDFs using Inverse-Gamma model (stars)
for r_idx = 1:length(r_vec)
    for k_idx = 1:length(K_u_vec)
        % Extract all SINRs
        SINRs = squeeze(SINR_RZF(r_idx, k_idx, :));
        SINRs_dB = 10 * log10(SINRs);

        % Estimate DSk^2 and IUSInk + noise as average (approximation)
        avg_SINR = mean(SINRs);
        avg_desired = mean(abs(v_RZF' * h_desired)^2);  % Approximate from last loop
        avg_known_noise = mean(sum(abs(v_RZF' * h_known).^2) + noise_power * norm(v_RZF)^2);
        
        % Estimate unknown interference samples from SINR rearrangement
        IUI_samples = (avg_desired ./ SINRs) - avg_known_noise;

        % Remove negative values (numerical noise)
        IUI_samples = IUI_samples(IUI_samples > 0);

        % Fit inverse-gamma to IUI_samples
        mu = mean(IUI_samples);
        v = var(IUI_samples);
        alpha = (mu^2 / v) + 2;
        beta = ((mu^2 / v) + 1) * mu;

        % Compute CDF analytically from inverse-gamma
        x_vals_dB = linspace(min(SINRs_dB), max(SINRs_dB), 300);
        x_vals = 10.^(x_vals_dB / 10);  % Convert dB to linear
        F_analytical = zeros(size(x_vals));
        for i = 1:length(x_vals)
            T = x_vals(i);
            gamma_input = (avg_desired / T) - avg_known_noise;
            if gamma_input > 0 
                % Inverse-Gamma CDF: F(x) = gammainc(beta / x, alpha)
                F_analytical(i) = gamcdf(beta / gamma_input, alpha,1);
            else
                F_analytical(i) = 0;
            end
        end

        % Overlay as stars
        plot(x_vals_dB, F_analytical, '*', ...
            'Color', colors(k_idx,:), ...
            'MarkerSize', 5, 'DisplayName', 'Analytical');
        %{
if r_idx == 1 && k_idx == 1
            PL = @(d) 10.^(-(30.5 + 36.7*log10(d))); % No shadowing for analytical
            % Average channel gains
            beta_desired = PL(r_vec(r_idx));            % Desired UE
            beta_known = PL([60, 100, 140, 180, 220]); % Known interferers
            beta_unknown = PL(250 + (500-250)*rand(1, 20)); % Unknown interferers
            
            % Total interference power
            I_total = P_tx * (sum(beta_known) + sum(beta_unknown)) + sigma2;
            
            % Mean and variance of SINR
            E_gamma = (P_tx * beta_desired) / I_total;
            Var_gamma = (P_tx^2 * beta_desired^2) / (N * I_total^2);
            
            % Gamma distribution parameters
            k = (E_gamma^2) / Var_gamma;     % Shape parameter
            theta = Var_gamma / E_gamma;     % Scale parameter
            
            % Outage probability (CDF of Gamma)
            T_dB = -30:0.1:20;               % SINR threshold range (dB)
            T_linear = 10.^(T_dB/10);        % Convert to linear
            P_out = gammainc(T_linear / theta, k); % CDF
            
            plot(T_dB, P_out, 'LineWidth', 2, 'DisplayName', 'Analytical MR');
        end
%}
    end
end
