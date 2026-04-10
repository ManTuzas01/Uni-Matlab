clc

%% INPUTS

m_cr = 2500; %[kg]
RPM = 395;
p = 27.8; %[kg/m^2]
C = 0.33; %[m]
m_ment = 11; %[kg]
n_ment = 3;
n_coef = 0.9;
v_max = 260; %[m/s]

g = 9.81; %[kg]

%% PRADINIAI SKAICIAVIMAI

% Maksimalus mentes ilgis

w = (2*pi*RPM)/60; %[rad/s]

R_max = v_max / w; %[m]

% Sraigto disko plotas

W_cr = m_cr * g; %[kg]

F = W_cr / (p * g); %[m^2]

% Mentes ilgis

R = sqrt(F/pi) %[m]

if R > R_max
    fprintf('Error R > R_max')
else
    fprintf('R < R_max == %4.2f m < %4.2f m \n', R, R_max)
end

% Oro srautas aplink mente

V = w * R; %[m/s]

fprintf('R_max = %4.2f m \n', R_max)
fprintf('F = %4.2f m \n', F)
fprintf('R = %4.2f m \n', R)
fprintf('V = %4.2f m \n', V)

%% KABOJIMUI REIKALINGA GALIA IR JEGA

% Traukos jega
T = W_cr; %[N]
T_ment = T/n_ment; %[N]

% Galia

V_vid = (w * R)/2; %[m/s]

P_ideal = (T * V_vid)/ 1000; %[kW]

P_prof = 0.25 * P_ideal; %[kW]

P = P_ideal + P_prof; %[kW]

P_vel = P/n_coef; %[kW]

P_gal = 1.2 * P_vel; %[kW]

fprintf('T = %4.2f N \n', T)
fprintf('T_ment = %4.2f N \n', T_ment)
fprintf('V_vid = %4.2f m/s \n', V_vid)
fprintf('P_ideal = %4.2f kW \n', P_ideal)
fprintf('P_prof = %4.2f kW \n', P_prof)
fprintf('P = %4.2f kW \n', P)
fprintf('P_vel = %4.2f kW \n', P_vel)
fprintf('P_gal = %4.2f kW \n', P_gal)

%% IRAZU SKAICIAVIMAS

%Saknyje
F_c = 0.5 * m_ment * w^2 * R;   % [N]

% Tiesinis mentes mases tankis
mu = m_ment / R; %[kg/m]

% Koordinate per mente nuo saknies iki galo
x = linspace(0, R, 500); %[m]

% Asine jega nuo iscentriniu jegu
Nx = 0.5 * mu * w^2 .* (R^2 - x.^2); %[N]

% Saknyje didziausia reiksme
Nx_root = Nx(1); %[N]

fprintf('mu = %4.4f kg/m \n', mu)
fprintf('Nx_root = %4.2f N \n', Nx_root)

figure;
plot(x, Nx, 'LineWidth', 2)
grid on
xlabel('Atstumas nuo mentes saknies x [m]')
ylabel('Asine jega N_x [N]')
title('Asines jegos pasiskirstymas menteje')

