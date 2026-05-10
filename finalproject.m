% CHEME 5440
% Final Project 

clear all; close all;
%% 

% in M 
% x0(1) = 16 * 10^-9; % initial concentration of free Ab (M), from figure on textbook, selected middle value
% x0(2) = 7.5 * 10^-11; % initial concentration of free CD4 receptor 
% x0(3) = 7.5 * 10^-11; % initial concentration of free CD70 receptor
% x0(4) = 0; % initial concentration of Ab:CD4
% x0(5) = 0; % initial concentration of Ab:CD70
% x0(6) = 0; % initial concentration of CD4:Ab:CD70
% 
% k(1) = 2.8 * 10^5; % kon, CD4 (M^-1 s^-1)
% k(2) = 2.6 * 10^-4; % koff, CD4 (s^-1)
% k(3) = 2.0 * 10^5; % kon, CD70 (M^-1 s^-1)
% k(4) = 4.9 * 10^-3; % koff, CD70 (s^-1)
% kappa = 2*10^6; % avidity parameter, from textbook 

% in nM 
% x0(1) = 16; % initial concentration of free Ab (nM)
% x0(2) = 7.5 * 10^-2; % initial conc of free CD4 receptor (nM)
% x0(3) = 7.5 * 10^-2; % initial conc of free CD70 receptor (nM)
% x0(4) = 0; % initial concentration of Ab:CD4 (nM)
% x0(5) = 0; % initial concentration of Ab:CD70 (nM)
% x0(6) = 0; % initial concentration of CD4:Ab:CD70 (nM)
% k(1) = 2.8 * 10^-4; % kon, CD4 (nM^-1 s^-1)
% k(2) = 2.6 * 10^-4; % koff, CD4 (s^-1)
% k(3) = 2.0 * 10^-4; % kon, CD70 (nM^-1 s^-1)
% k(4) = 4.9 * 10^-3; % koff, CD70 (s^-1)
% kappa = 2*10^6; % avidity parameter, from textbook

% Messing around in nM 
x0(1) = 16; % initial concentration of free Ab (nM)
x0(2) = 7.5 * 10^-2; % initial conc of free CD4 receptor (nM)
x0(3) = 7.5 * 10^-2; % initial conc of free CD70 receptor (nM)
x0(4) = 0; % initial concentration of Ab:CD4 (nM)
x0(5) = 0; % initial concentration of Ab:CD70 (nM)
x0(6) = 0; % initial concentration of CD4:Ab:CD70 (nM)
k(1) = 2.8 * 10^-4; % kon, CD4 (nM^-1 s^-1)
k(2) = 2.6 * 10^-4; % koff, CD4 (s^-1)
% Kd CD4 = around 0.9 nM 
k(3) = 2.0 * 10^-4; % kon, CD70 (nM^-1 s^-1)
k(4) = 4.9 * 10^-3; % koff, CD70 (s^-1)
% Kd CD70 = around 25 nM 
kappa = 2*10^2; % avidity parameter, from textbook

% Time-span
tspan = [0 60*60]; % time-span in sec

[t_out,x_out] = ode15s(@(t,x) finalprojectODEs(t,x,k,kappa),tspan,x0);

%Redfining CD70
k(1) = 5.2 * 10^-4; % kon, CD4 (nM^-1 s^-1)
k(2) = 3.6 * 10^-2; % koff, CD4 (s^-1)
% Kd CD70 = around 70 nM 

[t_out1,x_out1] = ode15s(@(t,x) finalprojectODEs(t,x,k,kappa),tspan,x0);

figure(1); 
subplot(2,3,1);
plot(t_out./60,x_out(:,1),'-r');
xlabel('Time (min)');
ylabel('Concentration (nM)');
title('antibody');
grid on;

subplot(2,3,2);
plot(t_out./60,x_out(:,2),'-r');
xlabel('Time (min)');
ylabel('Concentration (nM)');
title('cd4 receptor');
grid on;

subplot(2,3,3);
plot(t_out./60,x_out(:,3),'-r');
xlabel('Time (min)');
ylabel('Concentration (nM)');
title('cd70');
grid on;

subplot(2,3,4);
plot(t_out./60,x_out(:,4),'-r');
xlabel('Time (min)');
ylabel('Concentration (nM)');
title('Ab:CD4');
grid on;

subplot(2,3,5);
plot(t_out./60,x_out(:,5),'-r');
xlabel('Time (min)');
ylabel('Concentration (nM)');
title('Ab:CD70');
grid on;

subplot(2,3,6);
plot(t_out./60,x_out(:,6),'-r');
xlabel('Time (min)');
ylabel('Concentration (nM)');
title('CD4:Ab:CD70');
grid on;

figure(2); 
subplot(2,3,1);
plot(t_out./60,x_out(:,1),'-bo');
xlabel('Time (min)');
ylabel('Concentration (nM)');
title('antibody');
grid on;

subplot(2,3,2);
plot(t_out./60,x_out(:,2),'-bo');
xlabel('Time (min)');
ylabel('Concentration (nM)');
title('cd4 receptor');
grid on;

subplot(2,3,3);
plot(t_out./60,x_out(:,3),'-bo');
xlabel('Time (min)');
ylabel('Concentration (nM)');
title('cd70');
grid on;

subplot(2,3,4);
plot(t_out./60,x_out(:,4),'-bo');
xlabel('Time (min)');
ylabel('Concentration (nM)');
title('Ab:CD4');
grid on;

subplot(2,3,5);
plot(t_out./60,x_out(:,5),'-bo');
xlabel('Time (min)');
ylabel('Concentration (nM)');
title('Ab:CD70');
grid on;

subplot(2,3,6);
plot(t_out./60,x_out(:,6),'-bo');
xlabel('Time (min)');
ylabel('Concentration (nM)');
title('CD4:Ab:CD70');
grid on;

figure(3);
hold on;
plot(t_out./60,x_out(:,4),'-bo');
plot(t_out./60,x_out(:,5),'-go');
plot(t_out./60,x_out(:,6),'-ro');
xlabel('Time (min)');
ylabel('Concentration (nM)');
legend('Ab:CD4','Ab:CD70','CD4:Ab:CD70');
grid on;
hold off; 

figure(4);
for i = 1:6
    subplot(2,3,i);
    plot(t_out1./60,x_out1(:,i),'-r');
    xlabel('Time (min)');
    ylabel('Concentration (nM)');
    % title('CD4:Ab:CD70');
    grid on;
end

figure(5);
hold on;
plot(t_out1./60,x_out1(:,4),'-bo');
plot(t_out1./60,x_out1(:,5),'-go');
plot(t_out1./60,x_out1(:,6),'-ro');
xlabel('Time (min)');
ylabel('Concentration (nM)');
legend('Ab:CD4','Ab:CD70','CD4:Ab:CD70');
title(" now its more specific hopefully");
grid on;
hold off; 

% --- Setup Global Constants (Add these if not defined yet) ---
tspan = [0 5e5];        % Time for equilibrium
conv_factor = 6e5;      % Scale to #/cell
Ab_range = logspace(-11, -7, 25); 
kd_scenarios = [0.9, 70]; 
results = cell(1, 2);

% --- 1. Simulation Loop ---
for s = 1:2
    % Set kinetic rates
    k(1) = 2.8 * 10^-4;              % kon, CD4
    k(2) = k(1) * kd_scenarios(s);   % koff, CD4 
    k(3) = 2.0 * 10^-4;              % kon, CD70
    k(4) = k(3) * 25;                % koff, CD70 
    kappa = 2*10^2;                  
    
    bound_data = zeros(length(Ab_range), 3); 

    for i = 1:length(Ab_range)
        Ab_init_nM = Ab_range(i) * 1e9;
        
        % CD4+/CD70+ 
        x0 = [Ab_init_nM, 0.075, 0.075, 0, 0, 0];
        [~, x_out] = ode15s(@(t,x) finalprojectODEs(t,x,k,kappa), tspan, x0);
        bound_data(i, 1) = sum(x_out(end, 4:6)) * conv_factor;
        
        % CD4+/CD70-
        x0 = [Ab_init_nM, 0.075, 0, 0, 0, 0];
        [~, x_out] = ode15s(@(t,x) finalprojectODEs(t,x,k,kappa), tspan, x0);
        bound_data(i, 2) = sum(x_out(end, 4:6)) * conv_factor;
        
        % CD4-/CD70+
        x0 = [Ab_init_nM, 0, 0.075, 0, 0, 0];
        [~, x_out] = ode15s(@(t,x) finalprojectODEs(t,x,k,kappa), tspan, x0);
        bound_data(i, 3) = sum(x_out(end, 4:6)) * conv_factor;
    end
    results{s} = bound_data;
end

% --- 2. Plotting ---
figure(6); % Explicitly naming this Figure 6
set(gcf, 'Position', [100, 100, 1000, 400]); % Set window size

for s = 1:2
    subplot(1, 2, s);
    data = results{s};
    
    % Use semilogx for the [Ab] log scale
    semilogx(Ab_range, data(:,1), 'k-', 'LineWidth', 1.8); hold on;
    semilogx(Ab_range, data(:,2), 'k--', 'LineWidth', 1.2);
    semilogx(Ab_range, data(:,3), 'k:', 'LineWidth', 1.2);
    
    % Labels & Formatting
    xlabel('[Ab] (M)');
    ylabel('Bound Ab (#/cell)');
    % Fixed the num24 typo here:
    title(['K_{d,CD4} = ', num2str(kd_scenarios(s)), ' nM']);
    ylim([0 6e4]);
    grid on;
    
    % Dynamic labels placed near the end of the curves
    text(2e-8, data(end,1), 'CD4^+/CD70^+', 'FontSize', 8);
    text(2e-8, data(end,2), 'CD4^+/CD70^-', 'FontSize', 8);
    text(2e-8, data(end,3), 'CD4^-/CD70^+', 'FontSize', 8);
end