% CHEME 5440
% Final Project 

clear all; close all;
%% Recreating Textbook Figure on Relationship between Aeff and L

% Messing around in nM 
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


hold on
%For our updated scenario with a different [CD4eff] concentration
d_vals = linspace(10,60,6);
figure(1)
hold on
x = 0;
startlength = [1,4,9,14,19,25];
for i = 1:length(d_vals)
    d = d_vals(i);
    L_vals = linspace(startlength(i),80,1000); %Accounts for model breakdown after you get lower linke lengths
    Aeff_vals = zeros(size(L_vals)); 

    for j = 1:length(L_vals)
        L = L_vals(j);
    
        Aeff = 5.03 .* L.^(-3/2) .* exp((-6.58e-2 .* d.^2) ./ L) .* ( ...
            1 ...
            - (0.987 ./ L) ...
            + 0.139 .* (d.^2 ./ L.^2) ...
            - 2.51e-3 .* (d.^4 ./ L.^3) ...
            - (0.308 ./ L.^2) ...
            - 0.150 .* (d.^2 ./ L.^3) ...
            + 0.0204 .* (d.^4 ./ L.^4) ...
            - 5.17e-4 .* (d.^6 ./ L.^5) ...
            + 3.14e-6 .* (d.^8 ./ L.^6) ...
        );
    
        Aeff_vals(j) = Aeff*1e3;   % store scalar
    end
    semilogy(L_vals, Aeff_vals, 'LineWidth', 2);

end
%Updated model
set(gca, 'YScale', 'log');
ylim([1e-3 1e3]);
xlabel('Linker length L (aa)');
ylabel('[A_{eff}] (mM)');
title('Effective concentration vs linker length');
xlim([0 80])
legend('d=10 A', 'd = 20 A', 'd = 30 A',  'd = 40 A', 'd = 50 A', 'd = 60 A')
grid on;
hold off

%% Mazor et al. figure recreation

%Set up constants
tspan = [0 60*60];        % Time for equilibrium
conv_factor = 9e8;      % Scaling factor
Ab_range = logspace(-11, -7, 100); 
kd_scenarios = [0.9, 70]; %nM
results = cell(1, 2);
avagadro = 6.022e23;

%finding Aeff value
d = 22;
L = 15;
% calculating Aeff
Aeff = 5.03 * L^(-3/2) * exp((-6.58e-2 * d^2) / L) * ( ...
            1 ...
            - (0.987 / L) ...
            + 0.139 * (d.^2 / L^2) ...
            - 2.51e-3 * (d^4 / L.^3) ...
            - (0.308 / L^2) ...
            - 0.150 * (d^2 / L.^3) ...
            + 0.0204 * (d^4 / L^4) ...
            - 5.17e-4 * (d^6 / L^5) ...
            + 3.14e-6 * (d^8 / L^6) ...
)*1e9;

for s = 1:2
    % kinetic rates in molar
    k(1) = 2.8e-4; %converting nanomolars to molars         
    k(3) = 2.0e-4;
    k(2) = k(1) * kd_scenarios(s); 
    k(4) = k(3) * 25; 
    disp(k(2))
    
    bound_data = zeros(length(Ab_range), 3); 

    for i = 1:length(Ab_range)
        Ab_init_nM = Ab_range(i)*1e9;
        
        % CD4+/CD70+ 
        x0 = [Ab_init_nM, (4.6e4/avagadro)*conv_factor*1e9, (5.2e4/avagadro)*conv_factor*1e9, 0, 0, 0];
        [~, x_out] = ode15s(@(t,x) finalprojectODEsaltered(t,x,k,Aeff), tspan, x0);
        bound_data(i, 1) = (sum(x_out(end, 4:6)) * avagadro * 1e-9)/conv_factor;
            
        % CD4+/CD70-
        x0 = [Ab_init_nM, (3.8e4/avagadro)*conv_factor*1e9, 0, 0, 0, 0];
        [~, x_out] = ode15s(@(t,x) finalprojectODEsaltered(t,x,k,Aeff), tspan, x0);
        bound_data(i, 2) = (sum(x_out(end, 4:6)) * avagadro * 1e-9)/conv_factor;

        % CD4-/CD70+
        x0 = [Ab_init_nM, 0, (3.1e4/avagadro)*conv_factor*1e9, 0, 0, 0];
        [~, x_out] = ode15s(@(t,x) finalprojectODEsaltered(t,x,k,Aeff), tspan, x0);
        bound_data(i, 3) = (sum(x_out(end, 4:6)) * avagadro * 1e-9)/conv_factor;
    end
    
    results{s} = bound_data; 
end
set(gcf, 'Position', [100, 100, 1000, 400]); % Set window size

figure(2);
for s = 1:2
    subplot(1, 2, s);
    data = results{s};
    
    % Use semilogx for the [Ab] log scale
    semilogx(Ab_range, data(:,1), 'k-', 'LineWidth', 1.8); hold on;
    semilogx(Ab_range, data(:,2), 'k--', 'LineWidth', 1.2);
    semilogx(Ab_range, data(:,3), 'k:', 'LineWidth', 1.2);
    
    % Labels & Formatting
    xlabel('[scFvs] (M)');
    ylabel('Bound scFvs (#/cell)');
    % Fixed the num24 typo here:
    title(['K_{d,CD4} = ', num2str(kd_scenarios(s)), ' nM']);
    grid on;

    legend('CD4+/CD70+', 'CD4+/CD70-','CD4-/CD70+')

end

