clc
close all

% This script is intended to be run after running the Simulink model 'Experiments.slx'

% Set default figure and font settings
set(0, 'DefaultLineLineWidth', 2, 'DefaultAxesFontSize', 28, 'DefaultTextFontSize', 28);

% Plot true states and their estimates
figure;
subplot(2, 1, 1);
plot(t, w_m, t, w_m_hat_est, t,  w_m_hat_mid, '--', t,  w_m_hat_min, '-.');
% NOTE: The line of code above was corrected in a later version.
% Previously, w_m_hat_est was mistakenly plotted three times under different labels.
% The error only affects the state estimate subplots in the thesis as the performance indices and delay estimates were correct.
% The impact is purely visual and does not affect the conclusions.
grid on
ylabel('$\omega_m(t)$', 'Interpreter', 'latex');
legend('measurement', '$d^*(t) = \mathrm{sat}(\hat{d}(t^-))$', '$d^*(t) = (\underline{d}+\overline{d})/2$',...
'$d^*(t) = \underline{d}$', 'Interpreter', 'latex');

subplot(2, 1, 2);
p = plot(t, D1, t, D1_hat_est, t, D1_hat_mid, '--', t, D1_hat_min,'-.');
grid on
xlabel('time (s)', 'Interpreter', 'latex');
ylabel('$d(t)$', 'Interpreter', 'latex');
legend('real', '$d^*(t) = \mathrm{sat}(\hat{d}(t^-))$', '$d^*(t) = (\underline{d}+\overline{d})/2$',...
'$d^*(t) = \underline{d}$', 'Interpreter', 'latex');
ylim([0 1.25]); % Scenario 1,2,3,
%ylim([0 1.5]); % Scenario 4,6,9


% Plot derivative of input signals
default_colors = get(groot, 'defaultAxesColorOrder');
custom_colors = [default_colors(2,:); default_colors(3,:); default_colors(4:end,:)];
    
figure;
set(gcf, 'DefaultAxesColorOrder', custom_colors);
hold on;
plot(t, U1dot_est, t, U1dot_mid, '--', t, U1dot_min, '-.');
grid on
xlabel('time (s)', 'Interpreter', 'latex');
ylabel('$\dot{u}(t - d^{*}(t))$', 'Interpreter', 'latex');
legend('$d^*(t) = \mathrm{sat}(\hat{d}(t^-))$', '$d^*(t) = (\underline{d}+\overline{d})/2$',...
    '$d^*(t) = \underline{d}$', 'Interpreter', 'latex');


% Compute IAE and ISE of delays for t >= 1.5
start_idx = find(t >= 1.5, 1);
T_start = t(start_idx:end);

IAE_w_m_est = trapz(T_start, abs(w_m(start_idx:end) - w_m_hat_est(start_idx:end)));
ISE_w_m_est = trapz(T_start, (w_m(start_idx:end) - w_m_hat_est(start_idx:end)).^2);
IAE_d1_est = trapz(T_start, abs(D1(start_idx:end) - D1_hat_est(start_idx:end)));
ISE_d1_est = trapz(T_start, (D1(start_idx:end) - D1_hat_est(start_idx:end)).^2);

IAE_w_m_mid = trapz(T_start, abs(w_m(start_idx:end) - w_m_hat_mid(start_idx:end)));
ISE_w_m_mid = trapz(T_start, (w_m(start_idx:end) - w_m_hat_mid(start_idx:end)).^2);
IAE_d1_mid = trapz(T_start, abs(D1(start_idx:end) - D1_hat_mid(start_idx:end)));
ISE_d1_mid = trapz(T_start, (D1(start_idx:end) - D1_hat_mid(start_idx:end)).^2);

IAE_w_m_min = trapz(T_start, abs(w_m(start_idx:end) - w_m_hat_min(start_idx:end)));
ISE_w_m_min = trapz(T_start, (w_m(start_idx:end) - w_m_hat_min(start_idx:end)).^2);
IAE_d1_min = trapz(T_start, abs(D1(start_idx:end) - D1_hat_min(start_idx:end)));
ISE_d1_min = trapz(T_start, (D1(start_idx:end) - D1_hat_min(start_idx:end)).^2);
