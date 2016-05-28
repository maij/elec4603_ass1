%% Static parameters
N_A  = 1e17;                  % Acceptor concentration (cm^-3)
N_D  = 1e19;                  % Donor    concentration (cm^-3)
N_V  = 1.04e19;
N_C  = 2.8e19;          
n_i  = 1e10;                  % Intrinsic carrier concentration (cm^-3)

q    = 1.602e-19;             % Electron charge (magnitude)
k_B  = 1.381e-23;             % Boltzmann's constant (J / K )
T    = 300;                   % Temperature (K)

eps_0    = 8.854e-12;         % Electric permittivity of free space (F / m)
eps_si   = 11.68*eps_0;       % Electric permittivity of silicon


psi_0 = 1.12 - k_B*T*log(N_V*N_C/(N_A*N_D))/q
A_P = q*N_A/eps_si;
A_N = q*N_D/eps_si;

W = sqrt(2*eps_si/q*psi_0*(1/N_A + 1/N_D));

x_n = W/(1+N_D/N_A);
x_p = W - x_n;

x = linspace(-x_p, x_n, 1000);
% Calculate Psi
Psi = [];
for i = 1:1000
    if x(i) <= 0
        Psi(i) = q*N_A/eps_si * (0.5*x(i)^2 + x_p*x(i));
    else
        Psi(i) = - q*N_D/eps_si*(0.5*x(i)^2 - x_n*x(i));
    end
end
% Calculate E
E = [];
for i = 1:1000
    if x(i) <= 0
        E(i) = - q*N_A/eps_si * (x(i) + x_p);
    else
        E(i) = q*N_D/eps_si*(x(i) - x_n);
    end
end

%% Psi plotting
figure(1);
 plot(x/100*1e6,Psi+psi_0, 'LineWidth', 3)
% line([0 0],[0 min(Psi)], 'LineStyle', '-.', 'Color', 'k')
title('Built-in voltage of a pn diode');
xlabel('Distance from pn junction (\mum)');
ylabel('\psi (V)');
set(gca, 'FontSize', 18)
%% E plot
figure(2);
plot(x/100*1e6,E, 'LineWidth', 3)
line([0 0],[0 min(E)], 'LineStyle', '-.', 'Color', 'k')
title('Electric field across depletion region of pn diode');
xlabel('Distance from pn junction (\mum)');
ylabel('E (V/m)');
set(gca, 'FontSize', 18)
