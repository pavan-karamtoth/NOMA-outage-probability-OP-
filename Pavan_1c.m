  clc;
clear;
close all;

% System Configuration
N = 1e5;
snr_dB = 0:5:30;
snr_lin = 10.^(snr_dB / 10);

noisePower = 1;

% Power allocation
a_f = 0.85;               % Power to far user
a_n = 1 - a_f;        % Power to near user

% SNR decoding thresholds
threshold_far = 2;
threshold_near = 1;

% Uplink transmit powers
P_uf = 1;
P_un = 1;

% Channel variances
% For Far User
gain_SDf = 1;
gain_UfDf = 1;
gain_UnDf = 1;

% For Near User
gain_SDn = 1;
gain_UfDn = 1;
gain_UnDn = 1;

% Result arrays
outage_far_sim = zeros(size(snr_lin));
outage_far_theory = zeros(size(snr_lin));
outage_near_sim = zeros(size(snr_lin));
outage_near_theory = zeros(size(snr_lin));

% 2) Main SNR Loop 

for idx = 1:length(snr_lin)

    P_bs = snr_lin(idx);
    gamma_bs = P_bs / noisePower;
    gamma_Uf = P_uf / noisePower;
    gamma_Un = P_un / noisePower;

    % Far User
    h1 = sqrt(gain_SDf/2) * (randn(N,1) + 1j*randn(N,1));
    h2 = sqrt(gain_UfDf/2) * (randn(N,1) + 1j*randn(N,1));
    h3 = sqrt(gain_UnDf/2) * (randn(N,1) + 1j*randn(N,1));

    g1 = abs(h1).^2;
    g2 = abs(h2).^2;
    g3 = abs(h3).^2;

    SINR_far = (a_f * gamma_bs .* g1) ./ (a_n * gamma_bs .* g1 + gamma_Uf * g2 + gamma_Un * g3 + 1);
    outage_far_sim(idx) = mean(SINR_far < threshold_far);

    theta_f = threshold_far / ((a_f - threshold_far * a_n) * gamma_bs);
    t1 = gain_SDf + gain_UfDf * gamma_Uf * theta_f;
    t2 = gain_SDf + gain_UnDf * gamma_Un * theta_f;
    outage_far_theory(idx) = 1 - (gain_SDf^2 / (t1 * t2)) * exp(-theta_f / gain_SDf);

    % Near User 
    h4 = sqrt(gain_SDn/2) * (randn(N,1) + 1j*randn(N,1));
    h5 = sqrt(gain_UfDn/2) * (randn(N,1) + 1j*randn(N,1));
    h6 = sqrt(gain_UnDn/2) * (randn(N,1) + 1j*randn(N,1));

    g4 = abs(h4).^2;
    g5 = abs(h5).^2;
    g6 = abs(h6).^2;

    SINR_decode_xf = (a_f * gamma_bs .* g4) ./ (a_n * gamma_bs .* g4 + gamma_Uf * g5 + gamma_Un * g6 + 1);
    SINR_decode_xn = (a_n * gamma_bs .* g4) ./ (gamma_Uf * g5 + gamma_Un * g6 + 1);

    outage_near_sim(idx) = mean((SINR_decode_xf < threshold_far) | (SINR_decode_xn < threshold_near));

    theta1 = threshold_far / (gamma_bs * (a_f - threshold_far * a_n));
    theta2 = threshold_near / (gamma_bs * a_n);
    theta_n = max(theta1, theta2);

    t1_near = gain_SDn + gain_UfDn * gamma_Uf * theta_n;
    t2_near = gain_SDn + gain_UnDn * gamma_Un * theta_n;

    outage_near_theory(idx) = 1 - (gain_SDn^2 / (t1_near * t2_near)) * exp(-theta_n / gain_SDn);
end

% Plotting =
figure;

% Far user
semilogy(snr_dB, outage_far_theory, 'bs--', 'LineWidth', 1.5, 'MarkerSize', 4); hold on;
semilogy(snr_dB, outage_far_sim, 'bd-', 'LineWidth', 1.5, 'MarkerSize', 4);

% Near user
semilogy(snr_dB, outage_near_theory, 'md--', 'LineWidth', 1.5, 'MarkerSize', 4);
semilogy(snr_dB, outage_near_sim, 'go-', 'LineWidth', 1.5, 'MarkerSize', 4);

grid on;
xlabel('Transmit Power (dB)');
ylabel('Outage Probability');
title('Outage Probability of Near and Far Users in NOMA', 'FontSize', 14, 'FontWeight', 'bold');
legend({'Far User - Analytical', 'Far User - Simulation', ...
        'Near User - Analytical', 'Near User - Simulation'}, ...
        'Location', 'best', 'FontSize', 10);
set(gca, 'FontSize', 11);