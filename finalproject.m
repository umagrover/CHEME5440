% CHEME 5440
% Final Project 
cd('C:\Users\umaag\OneDrive\Cornell\CHEME 5440\')

clear all; close all;
%% 

% Using Units of nM 
x0(1) = 10; % initial concentration of free Ab (nM)
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
kappa = 2*10^6; % avidity parameter, from textbook

% Time-span
tspan = [0 60*60]; % time-span in sec

[t_out,x_out] = ode15s(@(t,x) finalprojectODEs(t,x,k,kappa),tspan,x0);
% Kd CD70 = around 70 nM

%% Recreating the Mazor et al. Figure
% Global Constants
tspan = [0 60*60];        % Time for equilibrium
conv_factor = 9e8; %300000/0.001; %cells/liter
Ab_range = logspace(-11, -7, 25); 
kd_scenarios = [0.9, 70]; 
results = cell(1, 2);
avagadro = 6.022e23;
kappa= 2e6;

for s = 1:2
    % Set kinetic rates
    k(1) = 2.8 * 10^-4;
    %k(1) = 2.8 * 10^-4;              % kon, CD4
    k(2) = k(1) * kd_scenarios(s);   % koff, CD4 
    k(3) = 2.0 * 10^-4;              % kon, CD70
    k(4) = k(3) * 25;                % koff, CD70 
    %kappa = (kd_scenarios(s)/25)*avagadro/conv_factor; %using the derivaition from the textbook

    bound_data = zeros(length(Ab_range), 3); 

    for i = 1:length(Ab_range)
        Ab_init_nM = Ab_range(i) * 1e9;
        
        % CD4+/CD70+ 
        x0 = [Ab_init_nM, (4.6e4/avagadro)*conv_factor*1e9, (5.2e4/avagadro)*conv_factor*1e9, 0, 0, 0];
        %x0 = [Ab_init_nM, 0.075, 0.075, 0, 0, 0];
        [~, x_out] = ode15s(@(t,x) finalprojectODEs(t,x,k,kappa), tspan, x0);
        bound_data(i, 1) = (sum(x_out(end, 4:6)) * avagadro * 1e-9)/conv_factor;
        
        % CD4+/CD70-
        x0 = [Ab_init_nM, (3.8e4/avagadro)*conv_factor*1e9, 0, 0, 0, 0];
        [~, x_out] = ode15s(@(t,x) finalprojectODEs(t,x,k,kappa), tspan, x0);
        bound_data(i, 2) = (sum(x_out(end, 4:6)) * avagadro * 1e-9)/conv_factor;
        
        % CD4-/CD70+
        x0 = [Ab_init_nM, 0, (3.1e4/avagadro)*conv_factor*1e9, 0, 0, 0];
        [~, x_out] = ode15s(@(t,x) finalprojectODEs(t,x,k,kappa), tspan, x0);
        bound_data(i, 3) = (sum(x_out(end, 4:6)) * avagadro * 1e-9)/conv_factor;
    end
    results{s} = bound_data;
end

% Plotting
figure(6); 
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
    
    legend('CD4+/CD70+', 'CD4+/CD70-','CD4-/CD70+')

end
