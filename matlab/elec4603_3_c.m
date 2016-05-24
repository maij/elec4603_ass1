%% Static parameters
N_A   = 1e19;                 % Acceptor donor concentration (cm^-2)
mu_n  = 1417;                 % Electron mobility (cm^2 / V sec)
tau_n = 25e-9;                % Carrier lifetime (s)
n_i   = 1e10;                 % Intrinsic carrier concentration (cm^-2)

n_p0 = 1e1;                   % Base minority carrier concentration, under T.E.
q    = 1.602e-19;             % Electron charge (magnitude)
k_B  = 1.381e-23;             % Boltzmann's constant (J / K )
T    = 300;                   % Temperature (K)

D_n  = mu_n * k_B * T / q;    % Diffusivity (cm^2 / sec)
L_n  = sqrt(D_n*tau_n);       % Minority carrier diffusion length (cm)

base_coeff     = n_i^2/N_A;   % Constant multiplier
exponent_coeff = q/(k_B * T); % Constant exponent

%% Test parameters
W_B  = 1e-3;                  % Width of base (cm)
V_BE = -0.1;                  % Voltage across base-emitter junction (V)
V_BC = -7;                    % Voltage across base-collector junction (V)

%% Simulation
x = linspace(0, W_B, 1e3);
delta_n_2 = base_coeff*(exp(exponent_coeff * (-V_BE))-1);
delta_n_3 = base_coeff*(exp(exponent_coeff * (V_BC))-1);

denom   = sinh(W_B/L_n);
% Calculate the minority carrier concentration for each x.
n_p = n_p0 + delta_n_2 * (sinh((W_B - x)/L_n)/denom) + delta_n_3 * (sinh(x/L_n)/denom);
low_recomb_approx  = n_p0 + delta_n_2 * (1 - x./W_B) + delta_n_3 *x./W_B;
high_recomb_approx = n_p0 + delta_n_2 * exp(- x./L_n) + delta_n_3 *exp((x-W_B)/L_n);

%% Plot results
figure(1);
semilogy(x*1e3, n_p, 'b-', 'LineWidth', 1.5);
set(gca, 'FontSize', 18);
line([min(x)*1e3 max(x)*1e3], [n_p0 n_p0], 'LineWidth', 1.5, 'LineStyle', '-.', 'Color', 'r');
title('Minority carrier concentration in base of npn transistor', ...
    'FontSize', 28);
xlabel('Distance from base-emitter junction (\mum)', 'FontSize', 28);
ylabel('Carrier concentration (cm^{-2})', 'FontSize', 28);
leg = legend('Minority carrier concentration', 'n_{p0}');
set(leg, 'FontSize', 18);

ylim_curr = get(gca,'ylim');
xlim([min(x), max(x)]*1e3);
ylim([min(1, ylim_curr(1)), ylim_curr(2)]);