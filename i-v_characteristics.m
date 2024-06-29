% Clear previous data
close all;
clear all;
clc;
% Given parameters
k_n = 120e-6; % A/V^2
V_tn = 0.6; % V
k_p = 100e-6; % A/V^2
V_tp = -0.6; % V
V_DD = 3.3; % V
% V_gs sweep values
V_gs_n = linspace(0, V_DD, 50);
V_gs_p = linspace(-V_DD, 0, 50);
% V_ds sweep values
V_ds = linspace(0, V_DD, 100);
V_ds_p = linspace(0, -V_DD, 100); % Define V_ds_p for pMOS

% Initialize current matrices
I_ds_n = zeros(length(V_gs_n), length(V_ds));
I_ds_p = zeros(length(V_gs_p), length(V_ds_p));

% Calculate nMOS I-V characteristics
for i = 1:length(V_gs_n)
    for j = 1:length(V_ds)
        if V_gs_n(i) < V_tn
            I_ds_n(i, j) = 0; % Cutoff region
        elseif V_ds(j) < (V_gs_n(i) - V_tn)
            I_ds_n(i, j) = k_n * ((V_gs_n(i) - V_tn) * V_ds(j) - 0.5 * V_ds(j)^2); % Triode region
        else
            I_ds_n(i, j) = 0.5 * k_n * (V_gs_n(i) - V_tn)^2; % Saturation region
        end
    end
end

% Calculate pMOS I-V characteristics
for i = 1:length(V_gs_p)
    for j = 1:length(V_ds_p)
        if V_gs_p(i) > V_tp
            I_ds_p(i, j) = 0; % Cutoff region
        elseif V_ds_p(j) > (V_gs_p(i) - V_tp)
            I_ds_p(i, j) = k_p * ((V_gs_p(i) - V_tp) * V_ds_p(j) - 0.5 * V_ds_p(j)^2); % Triode region
        else
            I_ds_p(i, j) = 0.5 * k_p * (V_gs_p(i) - V_tp)^2; % Saturation region
        end
    end
end

figure;
hold on;
for i = 1:5:length(V_gs_n)
    plot(V_ds, I_ds_n(i, :), 'DisplayName', sprintf('V_{gs} = %.2f V', V_gs_n(i)), 'LineWidth', 1.2);
end
grid on;
for i = 1:5:length(V_gs_p)
    plot(V_ds_p, -I_ds_p(i, :), 'DisplayName', sprintf('V_{gs} = %.2f V', V_gs_p(i)), 'LineWidth', 1.2);
end
xlabel('V (V)');
ylabel('I (A)');
title('I-V Characteristics');
grid on;
