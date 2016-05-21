%% Static parameters
N_A   = 1e19;                 % Acceptor donor concentration (cm^-2)
mu_n  = 1417;                 % Electron mobility (cm^2 / V sec)
tau_n = 25e-9;                % Carrier lifetime (s)
n_i   = 1e10;                % Intrinsic carrier concentration (cm^-2)

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
V_BE = 0.1;                   % Voltage across base-emitter junction (V)
V_CE = -7;                    % Voltage across base-collector junction (V)

%% Simulation
all_x = linspace(0, W_B, 50e3);
V_BC = V_BE - V_CE;
delta_n_2 = n_p0 + base_coeff*(exp(exponent_coeff * (-V_BE))-1);
delta_n_3 = n_p0 + base_coeff*(exp(exponent_coeff * (V_BC))-1);

common_denom   = sinh(W_B/L_n);

for i = 1:size(all_x,2)
    x = all_x(i);
    delta_n(i) = delta_n_2 * (sinh((W_B - x)/L_n)/common_denom) + delta_n_3 * (sinh(x/L_n)/common_denom);
end

%% Plot results
figure(1);
semilogy(all_x*1e3, delta_n, 'b-');
line([0 n_p0], [W_B n_p0], 'LineStyle', '-.', 'Color', 'r');
title('Minority carrier concentration in base of npn transistor');
xlabel('Distance from base-emitter junction (\mum)');
ylabel('Carrier concentration (cm^{-2})');
xlim([-1e-3, W_B*1e3]);
% plot(delta_n)