%% === Comprehensive Simulation: EHL vs Dry Contact (Inner-Race Spall) [SCI English Plotting] ===
clear; clc; close all;

global BPFI fr

%% 1. Unified simulation parameters
T_sim = 4;               % Simulation duration (s)
dt    = 1e-5;            % Sampling step (100 kHz)
tspan = 0:dt:T_sim;
y0    = [1e-6, 0, 1e-6, 0, 1e-6, 0, 1e-6, 0];
OPT   = odeset('RelTol', 1e-3, 'AbsTol', 1e-6);

% Frequency parameters
n = 1250; Dm = 72.5e-3; Db = 11.5e-3; Nb = 9;
fr = n / 60;
BPFI = (Nb / 2) * (1 + Db / Dm) * fr;

fprintf('Rotational Freq (fr) = %.2f Hz\n', fr);
fprintf('Fault Freq (BPFI)    = %.2f Hz\n', BPFI);

%% 2. Run both models
models = {'EHL', 'DRY'};
Data = struct();

for i = 1:length(models)
    mode = models{i};
    fprintf('\n>>> Running Model: %s ...\n', mode);

    tic;
    if strcmp(mode, 'EHL')
        [t, y] = ode45(@Bearing_Fault_EHL, tspan, y0, OPT);
    else
        [t, y] = ode45(@Bearing_Fault_Dry, tspan, y0, OPT);
    end
    sim_time = toc;
    fprintf('    Elapsed Time: %.2f s\n', sim_time);

    % --- Data post-processing ---
    % 1. Remove transient (first 0.1 s)
    cut_t = 0.1;
    idx = t > cut_t;
    t_cut = t(idx);
    y_cut = y(idx, :);

    % 2. Acceleration
    dt_vec = diff(t_cut);
    acc_y = diff(y_cut(:, 8)) ./ dt_vec;
    acc_y = [acc_y; acc_y(end)];

    % 3. Force reconstruction
    if strcmp(mode, 'EHL')
        [F_ball1, F_impact_total] = ReconstructForce_EHL(t_cut, y_cut);
    else
        [F_ball1, F_impact_total] = ReconstructForce_Dry(t_cut, y_cut);
    end

    % 4. Store data
    Data(i).name = mode;
    Data(i).t = t_cut;
    Data(i).acc = acc_y;
    Data(i).F_ball1 = F_ball1;
    Data(i).F_impact_total = F_impact_total;
    Data(i).y = y_cut;
end

%% =============================================================
%% 3. Plot comparison (SCI-style, English)
%% =============================================================

% Common font
std_font = 'Times New Roman';

% SCI colors
Color_EHL = [0, 0.4470, 0.7410]; % Blue
Color_Dry = [0.4940, 0.1840, 0.5560]; % Purple

% --- Fig 1: Time-domain acceleration ---
figure(1); set(gcf, 'Position', [0 420 400 300], 'Color', 'w', 'Name', 'Time Domain Comparison');

subplot(2, 1, 1);
plot(Data(1).t, Data(1).acc, 'Color', Color_EHL, 'LineWidth', 0.6);
title('(a) EHL Model: Acceleration Signal', 'FontName', std_font, 'FontSize', 11);
ylabel('Amp (m/s^2)', 'FontName', std_font);
xlim([Data(1).t(1), Data(1).t(1) + 0.2]); grid on;
set(gca, 'FontName', std_font);

subplot(2, 1, 2);
plot(Data(2).t, Data(2).acc, 'Color', Color_Dry, 'LineWidth', 0.6);
title('(b) Dry Contact Model: Acceleration Signal', 'FontName', std_font, 'FontSize', 11);
ylabel('Amp (m/s^2)', 'FontName', std_font);
xlim([Data(2).t(1), Data(2).t(1) + 0.2]); grid on;
xlabel('Time (s)', 'FontName', std_font);
set(gca, 'FontName', std_font);

% --- Fig 2: Amplitude detail ---
figure(2); set(gcf, 'Position', [400 420 400 300], 'Color', 'w', 'Name', 'Amplitude Detail');
hold on;
start_t = Data(1).t(1);
range_idx = Data(1).t < (start_t + 0.05);
plot(Data(1).t(range_idx), Data(1).acc(range_idx), 'Color', Color_EHL, 'LineWidth', 0.8);
plot(Data(2).t(range_idx), Data(2).acc(range_idx), 'Color', 'r', 'LineStyle', '--', 'LineWidth', 0.7);

legend({'EHL Model', 'Dry Contact Model'}, 'FontName', std_font);
title('Detailed Comparison of Acceleration Amplitude', 'FontName', std_font, 'FontSize', 11);
ylabel('Acceleration (m/s^2)', 'FontName', std_font);
xlabel('Time (s)', 'FontName', std_font);
grid on; box on; set(gca, 'FontName', std_font);

% --- Fig 3: Spectrum comparison ---
figure(3); set(gcf, 'Position', [700 300 500 450], 'Color', 'w', 'Name', 'Spectrum Comparison');
subplot(2, 1, 1);
[f1, p1] = FFT_Custom(Data(1).t, Data(1).acc);
plot(f1, p1, 'Color', Color_EHL);
title('(a) EHL Spectrum', 'FontName', std_font, 'FontSize', 11);
xlim([0, 3000]); grid on; ylabel('Amplitude', 'FontName', std_font);
set(gca, 'FontName', std_font);

subplot(2, 1, 2);
[f2, p2] = FFT_Custom(Data(2).t, Data(2).acc);
plot(f2, p2, 'Color', Color_Dry);
title('(b) Dry Spectrum', 'FontName', std_font, 'FontSize', 11);
xlim([0, 3000]); grid on;
ylabel('Amplitude', 'FontName', std_font); xlabel('Frequency (Hz)', 'FontName', std_font);
set(gca, 'FontName', std_font);

%% --- Fig 4: Envelope spectrum comparison ---
figure(4); set(gcf, 'Position', [0 0 600 500], 'Color', 'w', 'Name', 'Envelope Spectrum Comparison');
plot_colors = {Color_EHL, Color_Dry};
plot_titles = {'(a) EHL Model: Envelope Spectrum', '(b) Dry Contact Model: Envelope Spectrum'};

% Search half-window (Hz)
search_w = 3.0;

for i = 1:2
    subplot(2, 1, i);

    % 1. Envelope spectrum
    sig = Data(i).acc - mean(Data(i).acc);
    sig_env = abs(hilbert(sig));
    sig_env = sig_env - mean(sig_env);
    [fe, pe] = FFT_Custom(Data(i).t, sig_env);

    % 2. Plot
    plot(fe, pe, 'Color', plot_colors{i}, 'LineWidth', 1.0); hold on;

    % 3. Locate actual fi (BPFI)
    idx_range = find(fe > BPFI - 5 & fe < BPFI + 5);
    if ~isempty(idx_range)
        [~, max_idx] = max(pe(idx_range));
        real_fi = fe(idx_range(max_idx));
    else
        real_fi = BPFI;
    end

    % 4. Label targets [theory freq, priority, label]
    target_freqs = [ ...
        fr,              1, "f_r"; ...
        real_fi,         1, "f_i"; ...
        real_fi - fr,    0, "f_i - f_r"; ...
        real_fi + fr,    0, "f_i + f_r"; ...
        2 * real_fi,     1, "2f_i"; ...
        2 * real_fi - fr, 0, "2f_i - f_r"; ...
        2 * real_fi + fr, 0, "2f_i + f_r"; ...
        3 * real_fi,     1, "3f_i" ...
    ];

    % Y-axis max for arrow lengths
    y_limit_max = max(pe(fe > 10 & fe < 400));
    threshold_ratio = 0.05;

    % 5. Label peaks
    for k = 1:size(target_freqs, 1)
        f_theory = double(target_freqs(k, 1));
        label_str = target_freqs(k, 3);

        idx_local = find(fe > f_theory - search_w & fe < f_theory + search_w);

        if ~isempty(idx_local)
            [p_val, p_loc] = max(pe(idx_local));
            f_val = fe(idx_local(p_loc));

            if p_val > y_limit_max * threshold_ratio && f_val > 2
                if mod(k, 2) == 0
                    arrow_len = y_limit_max * 0.15;
                else
                    arrow_len = y_limit_max * 0.25;
                end

                start_y = p_val + arrow_len;
                end_y   = p_val;

                q = quiver(f_val, start_y, 0, end_y - start_y, 0, ...
                    'Color', 'k', 'LineWidth', 1.0, 'MaxHeadSize', 0.5);
                q.Clipping = 'off';

                text(f_val, start_y, ['$', char(label_str), '$'], ...
                    'Interpreter', 'latex', ...
                    'HorizontalAlignment', 'center', ...
                    'VerticalAlignment', 'bottom', ...
                    'FontSize', 11, ...
                    'FontName', std_font, ...
                    'Color', 'k');
            end
        end
    end

    % 6. Styling
    xlim([0, 350]);
    ylim([0, y_limit_max * 1.5]);

    title(plot_titles{i}, 'FontName', std_font, 'FontSize', 12);
    ylabel('Amplitude (m/s^2)', 'FontName', std_font);

    set(gca, 'FontName', std_font, 'TickDir', 'out', 'LineWidth', 1.0, 'Box', 'off');
    grid off;

    xl = xlim; yl = ylim;
    line(xl, [yl(2) yl(2)], 'Color', 'k', 'LineWidth', 1.0);
    line([xl(2) xl(2)], yl, 'Color', 'k', 'LineWidth', 1.0);
end
xlabel('Frequency (Hz)', 'FontName', std_font, 'FontSize', 11);

%% --- Print defect frequency table ---
fprintf('\n=======================================================\n');
fprintf('       Defect Frequencies & Modulation Sidebands       \n');
fprintf('=======================================================\n');
fprintf('Rotational Frequency (fr) = %.4f Hz\n', fr);
fprintf('Fault Frequency (BPFI)    = %.4f Hz\n', BPFI);
fprintf('-------------------------------------------------------\n');
fprintf('Harmonic    | Lower SB (-fr) | Center Freq | Upper SB (+fr)\n');
fprintf('-------------------------------------------------------\n');
fprintf('1x BPFI     | %8.4f Hz    | %8.4f Hz  | %8.4f Hz\n', BPFI - fr, BPFI, BPFI + fr);
fprintf('2x BPFI     | %8.4f Hz    | %8.4f Hz  | %8.4f Hz\n', 2 * BPFI - fr, 2 * BPFI, 2 * BPFI + fr);
fprintf('3x BPFI     | %8.4f Hz    | %8.4f Hz  | %8.4f Hz\n', 3 * BPFI - fr, 3 * BPFI, 3 * BPFI + fr);
fprintf('=======================================================\n');

%% --- Fig 999: Fig 19 reproduction ---
if exist('Data', 'var')
    [Fc_matrix, ~] = Reconstruct_Forces_For_Plots(Data(1).t, Data(1).y);
    t_force = Data(1).t;
    t_start_plot = 0.22; t_end_plot = 0.26;

    figure(999);
    set(gcf, 'Position', [150, 150, 800, 500], 'Color', 'w', 'Name', 'Fig19 Reproduction');

    hold on; box on;
    p1 = plot(t_force, Fc_matrix(:, 1), 'm-', 'LineWidth', 1);
    p2 = plot(t_force, Fc_matrix(:, 2), 'b-', 'LineWidth', 1);
    p3 = plot(t_force, Fc_matrix(:, 3), 'g-', 'LineWidth', 1);

    xlim([t_start_plot, t_end_plot]); ylim([0, max(max(Fc_matrix)) * 1.2]);
    xlabel('Time (s)', 'FontSize', 12, 'FontName', std_font);
    ylabel('Contact Force (N)', 'FontSize', 12, 'FontName', std_font);
    title('Dynamic Contact Forces Between Adjacent Rollers', 'FontSize', 13, 'FontName', std_font);
    legend([p1, p2, p3], {'Roller 1', 'Roller 2', 'Roller 3'}, 'Location', 'NorthEast', 'FontName', std_font);
    grid on; set(gca, 'FontName', std_font);

    idx_zoom = t_force >= t_start_plot & t_force <= t_end_plot;
    [max_val, max_idx_local] = max(Fc_matrix(idx_zoom, 2));
    global_indices = find(idx_zoom);
    t_peak = t_force(global_indices(max_idx_local));
    %#ok<NASGU>
    %#ok<NASGU>
    hold off;
end

disp('=== All plots generated successfully (English/SCI Style) ===');

%% =======================================================================
%% ========================= Core function definitions ====================
%% =======================================================================

function [f, P1] = FFT_Custom(t, x)
    L = length(t); Ts = mean(diff(t)); Fs = 1 / Ts;
    Y = fft(x); P2 = abs(Y / L);
    P1 = P2(1:floor(L / 2) + 1); P1(2:end-1) = 2 * P1(2:end-1);
    f = Fs * (0:(L / 2)) / L;
end

%% --- Force reconstruction: EHL ---
function [F_ball1_rec, F_impact_total] = ReconstructForce_EHL(t, y)
    % 1. Base parameters
    D = 90e-3; d = 55e-3; Dm = 72.5e-3; Db = 11.5e-3; Nb = 9;
    fi = 0.515; fo = 0.525; cr = 5e-6; E = 2.07e11; v = 0.3; mb = 0.002;
    n = 1250; w = n * pi / 30; L = 1.5e-3;

    % 2. Geometry and pre-computation
    Dpw = (D + d) / 2; gamma = Db / Dm; wc = 0.5 * (1 - gamma) * w;
    Di_race = Dm - Db; ro = Db / 2;

    % EHL and contact parameters
    alpha_oil = 2.3e-8; eta = 0.02;
    E1 = 2 / ((1 - v^2) / E + (1 - v^2) / E); % 2/E'

    % Curvature sums
    Ri = fi * Db; Ro = fo * Db;
    rou_in = 2 / Db + 2 / Db + 2 / (Dm - Db) - 1 / Ri;
    rou_out = 2 / Db + 2 / Db - 2 / (Dm + Db) - 1 / Ro;

    % --- Inner-ring Hertz ---
    Rxi = Db * (1 + gamma) / 2; Ryi = fi * Db / (2 * fi - 1); k_ratio_i = Ryi / Rxi;
    Ki1 = 1.0339 * (k_ratio_i)^0.636;
    F_ei = 1.5277 + 0.6013 * log(k_ratio_i);
    E_ei = 1.003 + 0.5968 * (1 / k_ratio_i);

    % --- Outer-ring Hertz ---
    Rxo = Db * (1 - gamma) / 2; Ryo = fo * Db / (2 * fo - 1); k_ratio_o = Ryo / Rxo;
    Ko1 = 1.0339 * (k_ratio_o)^0.636;
    F_eo = 1.5277 + 0.6013 * log(k_ratio_o);
    E_eo = 1.003 + 0.5968 * (1 / k_ratio_o);

    % --- Contact parameter a* ---
    a_star_i = (2 * Ki1^2 * E_ei / pi)^(1 / 3);
    a_star_o = (2 * Ko1^2 * E_eo / pi)^(1 / 3);

    % Stiffness coefficients
    e1 = sqrt(max(0, 1 - 1 / Ki1^2)); e2 = sqrt(max(0, 1 - 1 / Ko1^2));
    delta1i = 2 * F_ei / pi * ((1 - e1^2) * pi / (2 * E_ei))^(1 / 3);
    delta1o = 2 * F_eo / pi * ((1 - e2^2) * pi / 2 * E_eo)^(1 / 3);
    K_ci = (2 * sqrt(2) / 3) * (E / (1 - v^2)) * (delta1i)^(-3 / 2) * (rou_in)^(-1 / 2);
    K_ce = (2 * sqrt(2) / 3) * (E / (1 - v^2)) * (delta1o)^(-3 / 2) * (rou_out)^(-1 / 2);

    % Oil film stiffness reference
    C_r = 5000; Fr = 0.5 * C_r; Jr = 1 / Nb * (cosd(0)^2.5 + 2 * cosd(40)^2.5 + 2 * cosd(80)^2.5);
    Qmax = Fr * (cos(pi / 9))^1.5 / (Nb * Jr);
    U = (pi * (Dpw / Db / 2 * n * (1 - gamma^2)) * Db) / 60;
    U1 = U * eta / (E1 * Rxi); G = alpha_oil * E1;
    B1_base = 2.69 * U1^0.67 * G^0.53 * (Nb * Jr * E1)^0.067 * (cos(pi / Nb))^(-0.1005) * pi / Nb * (1 - 0.61 * exp(-0.73 * Ki1)) * Rxi^1.134;
    U2 = U * eta / (E1 * Rxo);
    B2_base = 2.69 * U2^0.67 * G^0.53 * (Nb * Jr * E1)^0.067 * cos(pi / Nb)^(-0.1005) * pi / Nb * (1 - 0.61 * exp(-0.73 * Ko1)) * Rxo^1.134;
    Koi_ref = (0.067 * B1_base * Fr^(-1.067))^(-1);
    Koe_ref = (0.067 * B2_base * Fr^(-1.067))^(-1);
    ho = 1e-6; hi = 1e-6;

    % 3. Time loop
    F_ball1_rec = zeros(length(t), 1); F_impact_total = zeros(length(t), 1);
    beta0 = pi / 2; phi_di = asin(L / Di_race); delta_max = (Db / 2) - sqrt((Db / 2)^2 - (L / 2)^2);
    Kb_nom = (1 / ((1 / K_ci)^(2 / 3) + (1 / K_ce)^(2 / 3)))^(3 / 2);
    term_is = max(0, 1 - (L / (2 * ro))^2); cos_theta_IS = sqrt(term_is);
    term_ie = max(0, 1 - (L / Di_race)^2); cos_theta_IE = sqrt(term_ie);
    denom_FIS = Dm * sqrt(4 * ro^2 - L^2);

    for k = 1:length(t)
        tk = t(k); yk = y(k, :); phi_is = mod(w * tk + beta0, 2 * pi);
        f_imp_sum = 0;

        for j = 1:Nb
            theta = wc * tk + 2 * pi * (j - 1) / Nb; phi_i = mod(theta, 2 * pi);
            angle_diff = phi_i - phi_is;
            while angle_diff > pi, angle_diff = angle_diff - 2 * pi; end
            while angle_diff < -pi, angle_diff = angle_diff + 2 * pi; end
            is_in = 0; h_in = 0;
            if abs(angle_diff) <= phi_di
                h_in = SpallDepth(angle_diff, phi_di, delta_max);
                if angle_diff > 0
                    is_in = 1;
                end
            end

            cd = 0.5 * cr * (1 - cos(abs(3 * pi / 2 - theta)));
            delta_j = (yk(1) - yk(5)) * cos(theta) + (yk(3) - yk(7)) * sin(theta) - cd - h_in - ho - hi;
            d_delta_j = (yk(2) - yk(6)) * cos(theta) + (yk(4) - yk(8)) * sin(theta);

            Fc = 0; F_imp = 0;

            if delta_j > 0
                Q_inst = Kb_nom * delta_j^(1.5); if Q_inst < 1e-5, Q_inst = 1e-5; end
                stiff_ratio = (Q_inst / Qmax)^0.067;
                K_oi_inst = Koi_ref * stiff_ratio; K_oe_inst = Koe_ref * stiff_ratio;
                K_i_inst = (1 / K_ci + 1 / K_oi_inst)^(-1); K_e_inst = (1 / K_ce + 1 / K_oe_inst)^(-1);
                Kt = (1 / ((1 / K_i_inst)^(2 / 3) + (1 / K_e_inst)^(2 / 3)))^(3 / 2);

                F_elastic = Kt * delta_j^(1.5);

                Q_real = max(F_elastic, 1e-3);
                a_i_dyn = 0.0236 * a_star_i * (Q_real / (rou_in))^(1 / 3);
                a_o_dyn = 0.0236 * a_star_o * (Q_real / (rou_out))^(1 / 3);

                C_i_dyn = 6 * pi * eta * Rxi^1.5 * a_i_dyn / sqrt(2) / hi^1.5;
                C_o_dyn = 6 * pi * eta * Rxo^1.5 * a_o_dyn / sqrt(2) / ho^1.5;

                C_eq = (C_i_dyn + C_o_dyn) / 4;

                F_sum = max(0, F_elastic + C_eq * d_delta_j);

                if is_in
                    term1 = (7 * w^2 * (Dm^2 - 4 * ro^2)^2) / (192 * Dm^2);
                    term2 = (5 / 6) * 9.8 * ro * (cos_theta_IE - cos_theta_IS);
                    term3 = (5 / 12) * Kt * (delta_j^2) / mb;
                    Vb2_sq = term1 + term2 + term3;
                    if Vb2_sq > 0
                        Vb2 = sqrt(Vb2_sq);
                        F_IS = (mb * w * (2 * Di_race * Dm - (Dm^2 - 4 * ro^2)) / denom_FIS) * Vb2;
                        F_imp = F_IS * cos_theta_IS;
                    end
                end
                Fc = F_sum;
            end

            f_imp_sum = f_imp_sum + F_imp;
            if j == 1, F_ball1_rec(k) = Fc; end
        end
        F_impact_total(k) = f_imp_sum;
    end
end

%% --- Force reconstruction: Dry ---
function [F_ball1_rec, F_impact_total] = ReconstructForce_Dry(t, y)
    D = 90e-3; d = 55e-3; Dm = 72.5e-3; Db = 11.5e-3; Nb = 9;
    fi = 0.515; fo = 0.525; cr = 5e-6; E = 2.07e11; v = 0.3; mb = 0.002;
    n = 1250; w = n * pi / 30; L = 1.5e-3;
    gamma = Db / Dm; wc = 0.5 * (1 - gamma) * w;
    Ri = fi * Db; Ro = fo * Db; rou_in = 2 / Db + 2 / Db + 2 / (Dm - Db) - 1 / Ri; rou_out = 2 / Db + 2 / Db - 2 / (Dm + Db) - 1 / Ro;
    delta_i = 0.577; delta_o = 0.68;
    kbi = (2 * sqrt(2) / 3) * (E / (1 - v^2)) * delta_i^(-1.5) * rou_in^(-0.5);
    kbo = (2 * sqrt(2) / 3) * (E / (1 - v^2)) * delta_o^(-1.5) * rou_out^(-0.5);
    kb = (1 / ((1 / kbi)^(2 / 3) + (1 / kbo)^(2 / 3)))^(1.5);

    beta0 = pi / 2; Di_race = Dm - Db; phi_di = asin(L / Di_race); delta_max = (Db / 2) - sqrt((Db / 2)^2 - (L / 2)^2);
    ro = Db / 2; term_is = max(0, 1 - (L / (2 * ro))^2); cos_theta_IS = sqrt(term_is);
    term_ie = max(0, 1 - (L / Di_race)^2); cos_theta_IE = sqrt(term_ie);
    denom_FIS = Dm * sqrt(4 * ro^2 - L^2);

    F_ball1_rec = zeros(length(t), 1); F_impact_total = zeros(length(t), 1);
    for k = 1:length(t)
        tk = t(k); yk = y(k, :); phi_is = mod(w * tk + beta0, 2 * pi); f_imp_sum = 0;
        for j = 1:Nb
            theta = wc * tk + 2 * pi * (j - 1) / Nb; phi_k = mod(theta, 2 * pi);
            angle_diff = phi_k - phi_is;
            while angle_diff > pi, angle_diff = angle_diff - 2 * pi; end
            while angle_diff < -pi, angle_diff = angle_diff + 2 * pi; end
            is_in = 0; h_in = 0;
            if abs(angle_diff) <= phi_di
                h_in = SpallDepth(angle_diff, phi_di, delta_max);
                if angle_diff > 0
                    is_in = 1;
                end
            end
            cd = 0.5 * cr * (1 - cos(abs(3 * pi / 2 - theta)));
            delta_j = (yk(1) - yk(5)) * cos(theta) + (yk(3) - yk(7)) * sin(theta) - cd - h_in;
            Fc = 0; F_imp = 0;
            if delta_j > 0
                F_elast = kb * delta_j^(1.5);
                if is_in
                    term1 = (7 * w^2 * (Dm^2 - 4 * ro^2)^2) / (192 * Dm^2);
                    term2 = (5 / 6) * 9.8 * ro * (cos_theta_IE - cos_theta_IS);
                    term3 = (5 / 12) * kb * (delta_j^2) / mb;
                    Vb2_sq = term1 + term2 + term3;
                    if Vb2_sq > 0
                        Vb2 = sqrt(Vb2_sq);
                        F_IS = (mb * w * (2 * Di_race * Dm - (Dm^2 - 4 * ro^2)) / denom_FIS) * Vb2;
                        F_imp = F_IS * cos_theta_IS;
                    end
                end
                Fc = F_elast;
            end
            f_imp_sum = f_imp_sum + F_imp;
            if j == 1, F_ball1_rec(k) = Fc; end
        end
        F_impact_total(k) = f_imp_sum;
    end
end

%% --- Model 1: EHL ---
function dy = Bearing_Fault_EHL(t, y)
    D = 90e-3; d = 55e-3; Dm = 72.5e-3; Db = 11.5e-3; Nb = 9;
    fi = 0.515; fo = 0.525; cr = 5e-6; E = 2.07e11; v = 0.3;
    mi = 1; mo = 1; mb = 0.002;
    ci = 1000; co = 1000; ki = 1.4e6; ko = 1.6e6;
    n = 1250; w = n * pi / 30;
    Wx = 500; Wy = 0; L = 1.5e-3; e = 0;
    gamma = Db / Dm; wc = 0.5 * (1 - gamma) * w;
    Di_race = Dm - Db; ro = Db / 2;
    alpha_oil = 2.3e-8; eta = 0.02; E1 = 2 / ((1 - v^2) / E + (1 - v^2) / E);
    Dpw = (D + d) / 2;
    nb = Dpw / Db / 2 * n * (1 - gamma^2); U = (pi * nb * Db) / 60;
    Rxi = Db * (1 + gamma) / 2; Ryi = fi * Db / (2 * fi - 1); k_ratio_i = Ryi / Rxi;
    Ki1 = 1.0339 * (k_ratio_i)^0.636;
    U1 = U * eta / (E1 * Rxi); G = alpha_oil * E1;
    Q_nom = 2500 / Nb * 4.37; Wi = Q_nom / (E1 * Rxi^2);
    Hi0 = 2.69 * U1^0.67 * G^0.53 * Wi^(-0.067) * (1 - 0.61 * exp(-0.73 * Ki1));
    hi = Hi0 * Rxi; ho = hi;
    Ri = fi * Db; Ro = fo * Db;
    rou_in = 2 / Db + 2 / Db + 2 / (Dm - Db) - 1 / Ri; rou_out = 2 / Db + 2 / Db - 2 / (Dm + Db) - 1 / Ro;
    delta_i = 0.577; delta_o = 0.68;
    kbi = (2 * sqrt(2) / 3) * (E / (1 - v^2)) * delta_i^(-1.5) * rou_in^(-0.5);
    kbo = (2 * sqrt(2) / 3) * (E / (1 - v^2)) * delta_o^(-1.5) * rou_out^(-0.5);
    Kb_nom = (1 / ((1 / kbi)^(2 / 3) + (1 / kbo)^(2 / 3)))^(1.5);
    C_eq = 500;
    Fx_elastic_damp = 0; Fy_elastic_damp = 0; F_ITV_x = 0; F_ITV_y = 0;
    beta0 = pi / 2; phi_is = mod(w * t + beta0, 2 * pi);
    phi_di = asin(L / Di_race);
    delta_max = (Db / 2) - sqrt((Db / 2)^2 - (L / 2)^2);
    theta0 = 0;
    term_is = max(0, 1 - (L / (2 * ro))^2); cos_theta_IS = sqrt(term_is); sin_theta_IS = sqrt(1 - term_is);
    term_ie = max(0, 1 - (L / Di_race)^2); cos_theta_IE = sqrt(term_ie);
    denom_FIS = Dm * sqrt(4 * ro^2 - L^2);

    for j = 1:Nb
        theta = wc * t + 2 * pi * (j - 1) / Nb + theta0; phi_k = mod(theta, 2 * pi);
        angle_diff = phi_k - phi_is;
        while angle_diff > pi, angle_diff = angle_diff - 2 * pi; end
        while angle_diff < -pi, angle_diff = angle_diff + 2 * pi; end
        is_in = 0; h_in = 0;
        if abs(angle_diff) <= phi_di
            h_in = SpallDepth(angle_diff, phi_di, delta_max);
            if angle_diff > 0
                is_in = 1;
            end
        end
        cd = 0.5 * cr * (1 - cos(abs(3 * pi / 2 - theta)));
        delta_j = (y(1) - y(5)) * cos(theta) + (y(3) - y(7)) * sin(theta) - cd - h_in - ho - hi;
        d_delta_j = (y(2) - y(6)) * cos(theta) + (y(4) - y(8)) * sin(theta);
        if delta_j > 0
            F_elast = Kb_nom * delta_j^(1.5);
            F_damp = C_eq * d_delta_j;
            F_sum = max(0, F_elast + F_damp);
            Fx_elastic_damp = Fx_elastic_damp + F_sum * cos(theta);
            Fy_elastic_damp = Fy_elastic_damp + F_sum * sin(theta);
            if is_in
                term1 = (7 * w^2 * (Dm^2 - 4 * ro^2)^2) / (192 * Dm^2);
                term2 = (5 / 6) * 9.8 * ro * (cos_theta_IE - cos_theta_IS);
                term3 = (5 / 12) * Kb_nom * (delta_j^2) / mb;
                Vb2_sq = term1 + term2 + term3;
                if Vb2_sq > 0
                    Vb2 = sqrt(Vb2_sq);
                    F_IS = (mb * w * (2 * Di_race * Dm - (Dm^2 - 4 * ro^2)) / denom_FIS) * Vb2;
                    F_ITV = F_IS;
                    term_x = cos_theta_IS * sin(phi_is) + sin_theta_IS * cos(phi_is);
                    term_y = cos_theta_IS * cos(phi_is) + sin_theta_IS * sin(phi_is);
                    F_ITV_x = F_ITV_x + F_ITV * term_x;
                    F_ITV_y = F_ITV_y + F_ITV * term_y;
                end
            end
        end
    end
    dy = zeros(8, 1);
    dy(1) = y(2); dy(2) = (Wx - ci * y(2) - ki * y(1) - F_ITV_x - Fx_elastic_damp + mi * e * w^2 * cos(w * t)) / mi;
    dy(3) = y(4); dy(4) = (Wy - ci * y(4) - ki * y(3) - F_ITV_y - Fy_elastic_damp + mi * e * w^2 * sin(w * t) - mi * 9.81) / mi;
    Fx_tot = F_ITV_x + Fx_elastic_damp; Fy_tot = F_ITV_y + Fy_elastic_damp;
    dy(5) = y(6); dy(6) = (Fx_tot - co * y(6) - ko * y(5)) / mo;
    dy(7) = y(8); dy(8) = (Fy_tot - mo * 9.81 - co * y(8) - ko * y(7)) / mo;
end

%% --- Model 2: Dry ---
function dy = Bearing_Fault_Dry(t, y)
    D = 90e-3; d = 55e-3; Dm = 72.5e-3; Db = 11.5e-3; Nb = 9;
    fi = 0.515; fo = 0.525; cr = 5e-6; E = 2.07e11; v = 0.3;
    mi = 1; mo = 1; mb = 0.002;
    ci = 1000; co = 1000; ki = 1.4e6; ko = 1.6e6;
    n = 1250; w = n * pi / 30;
    Wx = 500; Wy = 0; L = 1.5e-3; e = 0;
    gamma = Db / Dm; wc = 0.5 * (1 - gamma) * w;
    Di_race = Dm - Db; ro = Db / 2;
    Ri = fi * Db; Ro = fo * Db;
    rou_in = 2 / Db + 2 / Db + 2 / (Dm - Db) - 1 / Ri; rou_out = 2 / Db + 2 / Db - 2 / (Dm + Db) - 1 / Ro;
    delta_i = 0.577; delta_o = 0.68;
    kbi = (2 * sqrt(2) / 3) * (E / (1 - v^2)) * delta_i^(-1.5) * rou_in^(-0.5);
    kbo = (2 * sqrt(2) / 3) * (E / (1 - v^2)) * delta_o^(-1.5) * rou_out^(-0.5);
    Kb_nom = (1 / ((1 / kbi)^(2 / 3) + (1 / kbo)^(2 / 3)))^(1.5);
    Fx_elastic = 0; Fy_elastic = 0; F_ITV_x = 0; F_ITV_y = 0;
    beta0 = pi / 2; phi_is = mod(w * t + beta0, 2 * pi);
    phi_di = asin(L / Di_race);
    delta_max = (Db / 2) - sqrt((Db / 2)^2 - (L / 2)^2);
    theta0 = 0;
    term_is = max(0, 1 - (L / (2 * ro))^2); cos_theta_IS = sqrt(term_is); sin_theta_IS = sqrt(1 - term_is);
    term_ie = max(0, 1 - (L / Di_race)^2); cos_theta_IE = sqrt(term_ie);
    denom_FIS = Dm * sqrt(4 * ro^2 - L^2);

    for j = 1:Nb
        theta = wc * t + 2 * pi * (j - 1) / Nb + theta0; phi_k = mod(theta, 2 * pi);
        angle_diff = phi_k - phi_is;
        while angle_diff > pi, angle_diff = angle_diff - 2 * pi; end
        while angle_diff < -pi, angle_diff = angle_diff + 2 * pi; end
        is_in = 0; h_in = 0;
        if abs(angle_diff) <= phi_di
            h_in = SpallDepth(angle_diff, phi_di, delta_max);
            if angle_diff > 0
                is_in = 1;
            end
        end
        cd = 0.5 * cr * (1 - cos(abs(3 * pi / 2 - theta)));
        delta_j = (y(1) - y(5)) * cos(theta) + (y(3) - y(7)) * sin(theta) - cd - h_in;
        if delta_j > 0
            F_sum = Kb_nom * delta_j^(1.5);
            Fx_elastic = Fx_elastic + F_sum * cos(theta);
            Fy_elastic = Fy_elastic + F_sum * sin(theta);
            if is_in
                term1 = (7 * w^2 * (Dm^2 - 4 * ro^2)^2) / (192 * Dm^2);
                term2 = (5 / 6) * 9.8 * ro * (cos_theta_IE - cos_theta_IS);
                term3 = (5 / 12) * Kb_nom * (delta_j^2) / mb;
                Vb2_sq = term1 + term2 + term3;
                if Vb2_sq > 0
                    Vb2 = sqrt(Vb2_sq);
                    F_IS = (mb * w * (2 * Di_race * Dm - (Dm^2 - 4 * ro^2)) / denom_FIS) * Vb2;
                    F_ITV = F_IS;
                    term_x = cos_theta_IS * sin(phi_is) + sin_theta_IS * cos(phi_is);
                    term_y = cos_theta_IS * cos(phi_is) + sin_theta_IS * sin(phi_is);
                    F_ITV_x = F_ITV_x + F_ITV * term_x;
                    F_ITV_y = F_ITV_y + F_ITV * term_y;
                end
            end
        end
    end
    dy = zeros(8, 1);
    dy(1) = y(2); dy(2) = (Wx - ci * y(2) - ki * y(1) - F_ITV_x - Fx_elastic + mi * e * w^2 * cos(w * t)) / mi;
    dy(3) = y(4); dy(4) = (Wy - ci * y(4) - ki * y(3) - F_ITV_y - Fy_elastic + mi * e * w^2 * sin(w * t) - mi * 9.81) / mi;
    Fx_tot = F_ITV_x + Fx_elastic; Fy_tot = F_ITV_y + Fy_elastic;
    dy(5) = y(6); dy(6) = (Fx_tot - co * y(6) - ko * y(5)) / mo;
    dy(7) = y(8); dy(8) = (Fy_tot - mo * 9.81 - co * y(8) - ko * y(7)) / mo;
end

%% --- Helper: force decomposition ---
function [F_ball1, F_elast_b1, F_imp_b1] = Reconstruct_Dry_Components(t, y)
    D = 90e-3; d = 55e-3; Dm = 72.5e-3; Db = 11.5e-3; Nb = 9;
    fi = 0.515; fo = 0.525; cr = 5e-6; E = 2.07e11; v = 0.3; mb = 0.002;
    n = 1250; w = n * pi / 30; L = 1.5e-3;
    gamma = Db / Dm; wc = 0.5 * (1 - gamma) * w;
    Ri = fi * Db; Ro = fo * Db; rou_in = 2 / Db + 2 / Db + 2 / (Dm - Db) - 1 / Ri; rou_out = 2 / Db + 2 / Db - 2 / (Dm + Db) - 1 / Ro;
    delta_i = 0.577; delta_o = 0.68;
    kbi = (2 * sqrt(2) / 3) * (E / (1 - v^2)) * delta_i^(-1.5) * rou_in^(-0.5);
    kbo = (2 * sqrt(2) / 3) * (E / (1 - v^2)) * delta_o^(-1.5) * rou_out^(-0.5);
    kb = (1 / ((1 / kbi)^(2 / 3) + (1 / kbo)^(2 / 3)))^(1.5);
    beta0 = pi / 2; Di_race = Dm - Db; phi_di = asin(L / Di_race);
    delta_max = (Db / 2) - sqrt((Db / 2)^2 - (L / 2)^2);
    ro = Db / 2; term_is = max(0, 1 - (L / (2 * ro))^2); cos_theta_IS = sqrt(term_is);
    term_ie = max(0, 1 - (L / Di_race)^2); cos_theta_IE = sqrt(term_ie);
    denom_FIS = Dm * sqrt(4 * ro^2 - L^2);

    N = length(t); F_ball1 = zeros(N, 1); F_elast_b1 = zeros(N, 1); F_imp_b1 = zeros(N, 1);

    for k = 1:N
        tk = t(k); yk = y(k, :); phi_is = mod(w * tk + beta0, 2 * pi);
        j = 1; theta = wc * tk + 2 * pi * (j - 1) / Nb; phi_k = mod(theta, 2 * pi);
        angle_diff = phi_k - phi_is;
        while angle_diff > pi, angle_diff = angle_diff - 2 * pi; end
        while angle_diff < -pi, angle_diff = angle_diff + 2 * pi; end
        is_in = 0; h_in = 0;
        if abs(angle_diff) <= phi_di
            h_in = SpallDepth(angle_diff, phi_di, delta_max);
            if angle_diff > 0
                is_in = 1;
            end
        end
        cd = 0.5 * cr * (1 - cos(abs(3 * pi / 2 - theta)));
        delta_j = (yk(1) - yk(5)) * cos(theta) + (yk(3) - yk(7)) * sin(theta) - cd - h_in;
        f_el = 0; f_im = 0;
        if delta_j > 0
            f_el = kb * delta_j^(1.5);
            if is_in
                term1 = (7 * w^2 * (Dm^2 - 4 * ro^2)^2) / (192 * Dm^2);
                term2 = (5 / 6) * 9.8 * ro * (cos_theta_IE - cos_theta_IS);
                term3 = (5 / 12) * kb * (delta_j^2) / mb;
                Vb2_sq = term1 + term2 + term3;
                if Vb2_sq > 0
                    Vb2 = sqrt(Vb2_sq);
                    F_IS = (mb * w * (2 * Di_race * Dm - (Dm^2 - 4 * ro^2)) / denom_FIS) * Vb2;
                    f_im = F_IS * cos_theta_IS;
                end
            end
        end
        F_ball1(k) = f_el + f_im; F_elast_b1(k) = f_el; F_imp_b1(k) = f_im;
    end
end

function delta_ball1 = ComputeDelta_Roller1_Dry(t, y)
    Dm = 72.5e-3; Db = 11.5e-3; Nb = 9;
    fi = 0.515; fo = 0.525; cr = 5e-6;
    n = 1250; w = n * pi / 30; L = 1.5e-3;
    gamma = Db / Dm; wc = 0.5 * (1 - gamma) * w;
    Di_race = Dm - Db; ro = Db / 2;

    beta0 = pi / 2;
    phi_di = asin(L / Di_race);
    delta_max = (Db / 2) - sqrt((Db / 2)^2 - (L / 2)^2);

    delta_ball1 = zeros(length(t), 1);

    for k = 1:length(t)
        tk = t(k); yk = y(k, :);
        phi_is = mod(w * tk + beta0, 2 * pi);

        j = 1;
        theta = wc * tk + 2 * pi * (j - 1) / Nb;
        angle_diff = mod(theta - phi_is + pi, 2 * pi) - pi;

        h_in = 0;
        if abs(angle_diff) <= phi_di
            h_in = SpallDepth(angle_diff, phi_di, delta_max);
        end

        cd = 0.5 * cr * (1 - cos(abs(3 * pi / 2 - theta)));
        delta_j = (yk(1) - yk(5)) * cos(theta) + (yk(3) - yk(7)) * sin(theta) ...
                  - cd - h_in;

        delta_ball1(k) = delta_j;
    end
end

function h_in = SpallDepth(angle_diff, phi_di, delta_max)
    % Sloped entry into the spall, with a flat bottom afterwards.
    % ramp_in_ratio controls the fraction of the defect angle used for entry.
    ramp_in_ratio = 0.25;
    u = (angle_diff + phi_di) / (2 * phi_di);
    u = min(max(u, 0), 1);
    ramp_end = ramp_in_ratio;

    if u <= ramp_end
        h_in = delta_max * (u / ramp_end);
    else
        h_in = delta_max;
    end
end

%% --- Helper: multi-ball force reconstruction ---
function [Fc_all, phi_i_vec_all] = Reconstruct_Forces_For_Plots(t, y)
    D = 90e-3; d = 55e-3; Dm = 72.5e-3; Db = 11.5e-3; Nb = 9;
    fi = 0.515; fo = 0.525; cr = 5e-6; E = 2.07e11; v = 0.3; mb = 0.002;
    n = 1250; w = n * pi / 30; L = 1.5e-3;
    Dpw = (D + d) / 2; gamma = Db / Dm; wc = 0.5 * (1 - gamma) * w;
    Di_race = Dm - Db; ro = Db / 2;
    C_r = 5000; Fr = 0.1 * C_r; Jr = 1 / Nb * (cosd(0)^2.5 + 2 * cosd(40)^2.5 + 2 * cosd(80)^2.5); Qmax = Fr * (cos(pi / 9))^1.5 / (Nb * Jr);
    alpha_oil = 2.3e-8; eta = 0.02; E1 = 2 / ((1 - v^2) / E + (1 - v^2) / E);
    Rxi = Db * (1 + gamma) / 2; Ryi = fi * Db / (2 * fi - 1); k_ratio_i = Ryi / Rxi;
    Ki1 = 1.0339 * (k_ratio_i)^0.636; e1 = sqrt(max(0, 1 - 1 / Ki1^2));
    F_ei = 1.5277 + 0.6013 * log(k_ratio_i); E_ei = 1.003 + 0.5968 * (1 / k_ratio_i); delta1i = 2 * F_ei / pi * ((1 - e1^2) * pi / (2 * E_ei))^(1 / 3);
    U = (pi * (Dpw / Db / 2 * n * (1 - gamma^2)) * Db) / 60; U1 = U * eta / (E1 * Rxi); G = alpha_oil * E1;
    B1_base = 2.69 * U1^0.67 * G^0.53 * (Nb * Jr * E1)^0.067 * (cos(pi / Nb))^(-0.1005) * pi / Nb * (1 - 0.61 * exp(-0.73 * Ki1)) * Rxi^1.134;
    Rxo = Db * (1 - gamma) / 2; Ryo = fo * Db / (2 * fo - 1); k_ratio_o = Ryo / Rxo;
    Ko1 = 1.0339 * (k_ratio_o)^0.636; e2 = sqrt(max(0, 1 - 1 / Ko1^2));
    F_eo = 1.5277 + 0.6013 * log(k_ratio_o); E_eo = 1.003 + 0.5968 * (1 / k_ratio_o); delta1o = 2 * F_eo / pi * ((1 - e2^2) * pi / 2 * E_eo)^(1 / 3);
    U2 = U * eta / (E1 * Rxo); B2_base = 2.69 * U2^0.67 * G^0.53 * (Nb * Jr * E1)^0.067 * cos(pi / Nb)^(-0.1005) * pi / Nb * (1 - 0.61 * exp(-0.73 * Ko1)) * Rxo^1.134;
    Ri = fi * Db; Ro = fo * Db; rou_in = 2 / Db + 2 / Db + 2 / (Dm - Db) - 1 / Ri; rou_out = 2 / Db + 2 / Db - 2 / (Dm + Db) - 1 / Ro;
    K_ci = (2 * sqrt(2) / 3) * (E / (1 - v^2)) * (delta1i)^(-3 / 2) * (rou_in)^(-1 / 2);
    K_ce = (2 * sqrt(2) / 3) * (E / (1 - v^2)) * (delta1o)^(-3 / 2) * (rou_out)^(-1 / 2);
    Koi_ref = (0.067 * B1_base * Fr^(-1.067))^(-1); Koe_ref = (0.067 * B2_base * Fr^(-1.067))^(-1);
    ho = 1e-6; hi = 1e-6;
    roi_alli = 1 / Db * (4 - 1 / fi + 2 * gamma / (1 + gamma)); roi_allo = 1 / Db * (4 - 1 / fo - 2 * gamma / (1 - gamma)); Wx = 500;
    a1_damp = (2 * E_ei / (pi * (1 - e1^2)))^(1 / 3) * (3 * Wx / (2 * E1 * roi_alli))^(1 / 3);
    a_damp = (2 * E_eo / (pi * (1 - e1^2)))^(1 / 3) * (3 * Wx / (2 * E1 * roi_allo))^(1 / 3);
    Coi_base = 6 * pi * eta * Rxi^1.5 * a1_damp / sqrt(2) / hi^1.5; Coo_base = 6 * pi * eta * Rxo^1.5 * a_damp / sqrt(2) / ho^1.5;
    beta0 = pi / 2; phi_di = asin(L / Di_race); delta_max = (Db / 2) - sqrt((Db / 2)^2 - (L / 2)^2);
    Kb_nom = (1 / ((1 / K_ci)^(2 / 3) + (1 / K_ce)^(2 / 3)))^(3 / 2);
    term_is = max(0, 1 - (L / (2 * ro))^2); cos_theta_IS = sqrt(term_is);
    term_ie = max(0, 1 - (L / Di_race)^2); cos_theta_IE = sqrt(term_ie);
    denom_FIS = Dm * sqrt(4 * ro^2 - L^2);
    N_steps = length(t); Fc_all = zeros(N_steps, Nb); phi_i_vec_all = zeros(N_steps, Nb);
    for k = 1:N_steps
        tk = t(k); yk = y(k, :); phi_is = mod(w * tk + beta0, 2 * pi);
        for j = 1:Nb
            theta = wc * tk + 2 * pi * (j - 1) / Nb; phi_i = mod(theta, 2 * pi); phi_i_vec_all(k, j) = phi_i;
            angle_diff = phi_i - phi_is;
            while angle_diff > pi, angle_diff = angle_diff - 2 * pi; end
            while angle_diff < -pi, angle_diff = angle_diff + 2 * pi; end
            is_in = 0; h_in = 0;
            if abs(angle_diff) <= phi_di
                h_in = SpallDepth(angle_diff, phi_di, delta_max);
                if angle_diff > 0
                    is_in = 1;
                end
            end
            cd = 0.5 * cr * (1 - cos(abs(3 * pi / 2 - theta)));
            delta_j = (yk(1) - yk(5)) * cos(theta) + (yk(3) - yk(7)) * sin(theta) - cd - h_in - ho - hi;
            d_delta_j = (yk(2) - yk(6)) * cos(theta) + (yk(4) - yk(8)) * sin(theta);
            Fc_total = 0;
            if delta_j > 0
                Q_inst = Kb_nom * delta_j^(1.5); if Q_inst < 1e-5, Q_inst = 1e-5; end
                stiff_ratio = (Q_inst / Qmax)^0.067;
                K_oi_inst = Koi_ref * stiff_ratio; K_oe_inst = Koe_ref * stiff_ratio;
                K_i_inst = (1 / K_ci + 1 / K_oi_inst)^(-1); K_e_inst = (1 / K_ce + 1 / K_oe_inst)^(-1);
                Kt = (1 / ((1 / K_i_inst)^(2 / 3) + (1 / K_e_inst)^(2 / 3)))^(3 / 2);
                F_elastic = Kt * delta_j^(1.5);
                C_eq = (Coi_base + Coo_base) / 4;
                F_sum = max(0, F_elastic + C_eq * d_delta_j);
                F_imp_val = 0;
                if is_in
                    term1 = (7 * w^2 * (Dm^2 - 4 * ro^2)^2) / (192 * Dm^2);
                    term2 = (5 / 6) * 9.8 * ro * (cos_theta_IE - cos_theta_IS);
                    term3 = (5 / 12) * Kt * (delta_j^2) / mb;
                    Vb2_sq = term1 + term2 + term3;
                    if Vb2_sq > 0
                        Vb2 = sqrt(Vb2_sq);
                        F_IS = (mb * w * (2 * Di_race * Dm - (Dm^2 - 4 * ro^2)) / denom_FIS) * Vb2;
                        F_imp_val = F_IS * cos_theta_IS;
                    end
                end
                Fc_total = F_sum + F_imp_val;
            end
            Fc_all(k, j) = Fc_total;
        end
    end
end

%% === Fig 10: Global view ===
figure(10);
set(gcf, 'Position', [50, 100, 1200, 350], 'Color', 'w', 'Name', 'Global View (Ultra-Wide)');

Color_EHL = [0, 0.4470, 0.7410];

hold on;

p1 = plot(Data(1).t, Data(1).acc, 'Color', Color_EHL, 'LineWidth', 0.5);

p2 = plot(Data(2).t, Data(2).acc, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 0.5);

zoom_start = Data(1).t(1);
zoom_dur   = 0.045;
zoom_end   = zoom_start + zoom_dur;
%#ok<NASGU>

% Get Y-axis limits and set box height
%#ok<NASGU>
%#ok<NASGU>

xlim([Data(1).t(1), Data(1).t(end)]);

legend([p1, p2], {'EHL Model', 'Dry Contact Model'}, ...
    'FontName', 'Times New Roman', 'FontSize', 10, 'Location', 'northeast');

title(' Acceleration Signal Comparison ', 'FontName', 'Times New Roman', 'FontSize', 12);
ylabel('Acceleration (m/s^2)', 'FontName', 'Times New Roman');
xlabel('Time (s)', 'FontName', 'Times New Roman');

set(gca, 'FontName', 'Times New Roman', 'FontSize', 10, 'LineWidth', 1.0, 'TickDir', 'out');

xticks(0.1 : 0.05 : 0.4);

grid on;

hold off;

% --- Fig 5: Force vs Time ---
figure(5); set(gcf, 'Position', [0 420 400 300], 'Color', 'w', 'Name', 'Force vs Time Comparison');

subplot(2, 1, 1);
plot(Data(1).t, Data(1).F_ball1, 'Color', Color_EHL, 'LineWidth', 0.6);
title('(a) EHL Model: Force vs Time', 'FontName', std_font, 'FontSize', 11);
ylabel('Force (N)', 'FontName', std_font);
xlim([Data(1).t(1), Data(1).t(1) + 0.2]); grid on;
set(gca, 'FontName', std_font);

subplot(2, 1, 2);
plot(Data(2).t, Data(2).F_ball1, 'Color', Color_Dry, 'LineWidth', 0.6);
title('(b) Dry Contact Model: Force vs Time', 'FontName', std_font, 'FontSize', 11);
ylabel('Force (N)', 'FontName', std_font);
xlim([Data(2).t(1), Data(2).t(1) + 0.2]); grid on;
xlabel('Time (s)', 'FontName', std_font);
set(gca, 'FontName', std_font);

% --- Fig 6: Impact force vs Time ---
figure(6); set(gcf, 'Position', [400 420 400 300], 'Color', 'w', 'Name', 'Impact Force vs Time Comparison');

subplot(2, 1, 1);
plot(Data(1).t, Data(1).F_impact_total, 'Color', Color_EHL, 'LineWidth', 0.6);
title('(a) EHL Model: Impact Force vs Time', 'FontName', std_font, 'FontSize', 11);
ylabel('Impact Force (N)', 'FontName', std_font);
xlim([Data(1).t(1), Data(1).t(1) + 0.2]); grid on;
set(gca, 'FontName', std_font);

subplot(2, 1, 2);
plot(Data(2).t, Data(2).F_impact_total, 'Color', Color_Dry, 'LineWidth', 0.6);
title('(b) Dry Contact Model: Impact Force vs Time', 'FontName', std_font, 'FontSize', 11);
ylabel('Impact Force (N)', 'FontName', std_font);
xlim([Data(2).t(1), Data(2).t(1) + 0.2]); grid on;
xlabel('Time (s)', 'FontName', std_font);
set(gca, 'FontName', std_font);

disp('=== Force and Impact Force vs Time plots generated successfully ===');
