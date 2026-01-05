clc; 
clear;
close all;

% Simulation Configuration
N = 1e5;                           % Monte Carlo trials
snr_dB = 0:5:30;                   % BS transmit power in dB
snr_lin = 10.^(snr_dB / 10);       % Convert to linear scale

% System Parameters
noisePower = 1;
uplinkPower_far  = 5;
uplinkPower_near = 10;
powerAlloc_far   = 0.8;
powerAlloc_near  = 0.2;
threshold_far    = 2;
threshold_near   = 2;

% Channel Gains (mean of exponential RVs)
avgGain_BS_near   = 1;
avgGain_Uf_near   = 1;
avgGain_Un_near   = 1;

% Result Storage
outage_theory = zeros(size(snr_lin));
outage_sim    = zeros(size(snr_lin));

% Main Loop: Simulation & Theory
for idx = 1:length(snr_lin)
    
    % Power levels
    P_tx = snr_lin(idx);                    % Current BS power
    gamma_bs = P_tx / noisePower;           % SNR from BS
    gamma_uf = uplinkPower_far / noisePower;
    gamma_un = uplinkPower_near / noisePower;

    % (A) Analytical Outage 
    denom = powerAlloc_far - threshold_far * powerAlloc_near;

    if denom <= 0
        outage_theory(idx) = 1;
    else
        theta1 = threshold_far / (denom * gamma_bs);
        theta2 = threshold_near / (powerAlloc_near * gamma_bs);
        theta  = max(theta1, theta2);

        exp_term = exp(-theta / avgGain_BS_near);
        den1 = avgGain_BS_near + theta * gamma_uf * avgGain_Uf_near;
        den2 = avgGain_BS_near + theta * gamma_un * avgGain_Un_near;

        p_success = exp_term * (avgGain_BS_near^2 / (den1 * den2));
        outage_theory(idx) = 1 - p_success;
    end

    %  (B) Monte-Carlo Simulation
    % Generate complex Rayleigh fading samples

    h1 = sqrt(avgGain_BS_near/2) * (randn(N,1) + 1i*randn(N,1));
    h2 = sqrt(avgGain_Uf_near/2) * (randn(N,1) + 1i*randn(N,1));
    h3 = sqrt(avgGain_Un_near/2) * (randn(N,1) + 1i*randn(N,1));

    % Channel magnitudes (exponential)
    g1 = abs(h1).^2;
    g2 = abs(h2).^2;
    g3 = abs(h3).^2;

    % Compute SINRs
    SINR_decode_xf = (powerAlloc_far * P_tx .* g1) ./ ...
                     (powerAlloc_near * P_tx .* g1 + uplinkPower_far * g2 + uplinkPower_near * g3 + noisePower);
                 
    SINR_decode_xn = (powerAlloc_near * P_tx .* g1) ./ ...
                     (uplinkPower_far * g2 + uplinkPower_near * g3 + noisePower);

    % Outage event if either fails
    outages = (SINR_decode_xf < threshold_far) | (SINR_decode_xn < threshold_near);
    outage_sim(idx) = mean(outages);
end

% Plotting Results 

figure;
semilogy(snr_dB, outage_theory, 'md--', 'LineWidth', 1.5, 'MarkerSize', 4); hold on;
semilogy(snr_dB, outage_sim, 'go-', 'LineWidth', 1.5, 'MarkerSize', 4);
grid on;

xlabel('Transmit Power (dB)');
ylabel('Outage Probability');
title('Outage Probability at Near User (NOMA)', 'FontSize', 14, 'FontWeight', 'bold');
legend('Analytical', 'Monte Carlo', 'Location', 'best');
set(gca, 'FontSize', 11);
