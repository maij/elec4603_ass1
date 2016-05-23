%% Static parameters
N_A  = 1e19;                  % Acceptor donor concentration (cm^-3)
n_i  = 1e10;                  % Intrinsic carrier concentration (cm^-3)

N_A  = N_A * 100^3;           % Convert to (m^-3)
n_i  = n_i * 100^3;           % Convert to (m^-3)

q    = 1.602e-19;             % Electron charge (magnitude)
k_B  = 1.381e-23;             % Boltzmann's constant (J / K )
T    = 300;                   % Temperature (K)

eps_0    = 8.854e-12;         % Electric permittivity of free space (F / m)
eps_si   = 11.68*eps_0;       % Electric permittivity of silicon
eps_sio2 = 3.9*eps_0;         % Electric permittivity of silicon dioxide
d        = 1e-9;              % Dielectric width (m)
%% Stuff
phi_F = k_B*T/q*log(N_A/n_i)
Q_s = 2*sqrt(eps_si*q*N_A*phi_F)
Q_i = 0
C_i = eps_sio2/d

V_T = 2*phi_F - (Q_s+Q_i)/C_i