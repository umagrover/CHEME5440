% CHEME 5440
% Final Project 

clear all; close all;
%% Recreating Textbook Figure on Relationship between Aeff and L

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


hold on
%For our updated scenario iwth a different [CD4eff] concentration
d= 3.8;

d_vals = linspace(10,60,6);
figure (1)
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
conv_factor = 9e8;      % Scale to #/cell
Ab_range = logspace(-11, -2, 100); 
kd_scenarios = [0.9, 70]; 
results = cell(1, 2);
avagardo = 6.022e23;

%finding Aeff value
d = 40;
L = 40;
Aeff = (5.03 .* L.^(-3/2) .* exp((-6.58e-2 .* d.^2) ./ L) .* ( ...
            1 ...
            - (0.987 ./ L) ...
            + 0.139 .* (d.^2 ./ L.^2) ...
            - 2.51e-3 .* (d.^4 ./ L.^3) ...
            - (0.308 ./ L.^2) ...
            - 0.150 .* (d.^2 ./ L.^3) ...
            + 0.0204 .* (d.^4 ./ L.^4) ...
            - 5.17e-4 .* (d.^6 ./ L.^5) ...
            + 3.14e-6 .* (d.^8 ./ L.^6) ...
))*1e9;

%Loop
d = d_vals(4); %Using d = 40 A
for s = 1:2
    % Set kinetic rates
    k(1) = 2.8 * 10^-4;              % kon, CD4
    k(2) = k(1) * kd_scenarios(s);   % koff, CD4 
    k(3) = 2.0 * 10^-4;              % kon, CD70
    k(4) = k(3) * 25;                % koff, CD70            
    
    bound_data = zeros(length(Ab_range), 3); 

    pcd4 = 1;
    pcd70 =1;

    for i = 1:length(Ab_range)
        Ab_init_nM = Ab_range(i) * 1e9;
        
        o=50;
        % CD4+/CD70+ 
        x0 = [Ab_init_nM, 0.075, 0.075, 0, 0, 0];
        [~, x_out] = ode15s(@(t,x) finalprojectODEsaltered(t,x,k,Aeff_vals(o),pcd4,pcd70), tspan, x0);
        bound_data(i, 1) = (sum(x_out(end, 4:6)) * avagardo * 1e-9)/conv_factor;
            
        % CD4+/CD70-
        x0 = [Ab_init_nM, 0.075, 0, 0, 0, 0];
        [~, x_out] = ode15s(@(t,x) finalprojectODEsaltered(t,x,k,Aeff_vals(o),pcd4,pcd70), tspan, x0);
        bound_data(i, 2) = (sum(x_out(end, 4:6)) * avagardo * 1e-9)/conv_factor;
            
        % CD4-/CD70+
        x0 = [Ab_init_nM, 0, 0.075, 0, 0, 0];
        [~, x_out] = ode15s(@(t,x) finalprojectODEsaltered(t,x,k,Aeff_vals(o),pcd4,pcd70), tspan, x0);
        bound_data(i, 3) = (sum(x_out(end, 4:6)) * avagardo * 1e-9)/conv_factor;
    end
    results{s} = bound_data; 
end
xlim([1e-12,1e-6])
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
    xlabel('[Ab] (M)');
    ylabel('Bound Ab (#/cell)');
    % Fixed the num24 typo here:
    title(['K_{d,CD4} = ', num2str(kd_scenarios(s)), ' nM']);
    %ylim([0 6e4]);
    %xlim([10e-12,10e-8])
    grid on;
    
    % Dynamic labels placed near the end of the curves
    text(2e-8, data(end,1), 'CD4^+/CD70^+', 'FontSize', 8);
    text(2e-8, data(end,2), 'CD4^+/CD70^-', 'FontSize', 8);
    text(2e-8, data(end,3), 'CD4^-/CD70^+', 'FontSize', 8);
end
xlim([1e-12,1e-6])


%% Linker Length Variation and Monovalent/Bivalent Binding
%Bound complex 2 while varying linker length
Ab_vals_2 = [Ab_range(40), Ab_range(50), Ab_range(60)] ; %Choosing values to test
c2 = zeros(length(L_vals),length(Ab_vals_2));
c1 = zeros(length(L_vals),length(Ab_vals_2));
for c = 1:length(L_vals)
    L = L_vals(c);
    %Aeff in molars, so convert to nM, again using d = 40 from above
    Aeff = (5.03 .* L.^(-3/2) .* exp((-6.58e-2 .* d.^2) ./ L) .* ( ...
        1 ...
        - (0.987 ./ L) ...
        + 0.139 .* (d.^2 ./ L.^2) ...
        - 2.51e-3 .* (d.^4 ./ L.^3) ...
        - (0.308 ./ L.^2) ...
        - 0.150 .* (d.^2 ./ L.^3) ...
        + 0.0204 .* (d.^4 ./ L.^4) ...
        - 5.17e-4 .* (d.^6 ./ L.^5) ...
        + 3.14e-6 .* (d.^8 ./ L.^6) ...
     ))*1e9;

    for i = 1:length(Ab_vals_2)
        x0 = [Ab_vals_2(i), 0.075, 0.075, 0, 0, 0];
        [~, x_out2] = ode15s(@(t,x) finalprojectODEsaltered(t,x,k,Aeff,pcd4,pcd70), tspan, x0);
        c2(c,i) = (x_out2(end, 6) * avagardo * 1e-9)/conv_factor;
        c1(c,i) = (sum(x_out2(end, 4:5)) * avagardo * 1e-9)/conv_factor;
    end
end
figure(3)
subplot(1,2,1)
hold on
for o = 1:length(Ab_vals_2)
    plot(L_vals,c2(:,o),'Linewidth',2)

    %labels
    text(L_vals(end)*0.95, c2(end,o), num2str(Ab_vals_2(o)), 'FontSize', 8);
end
hold off
set(gca, 'YScale', 'log');
xlabel('Linker Length')
ylabel('C2 complex')
grid on
hold off

%c1
subplot(1,2,2)
hold on
for o = 1:length(Ab_vals_2)
    plot(L_vals,c1(:,o),'Linewidth',2)

    %labels
    text(L_vals(end)*0.95, c1(end,o), num2str(Ab_vals_2(o)), 'FontSize', 8);
end
hold off
set(gca, 'YScale', 'log');
xlabel('Linker Length')
ylabel('C1 complex')

grid on
hold off
