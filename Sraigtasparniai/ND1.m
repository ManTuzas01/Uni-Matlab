clc

%% INPUTS

m_cr = 2500; %[kg]
RPM = 395;
p = 27.8; %[kg/m^2]
C = 0.33; %[m]
m_ment = 11; %[kg/m]
n_ment = 3;
n_coef = 0.9;
v_max = 260; %[m/s]

g = 9.81; %[m/s^2]

%% PRADINIAI SKAICIAVIMAI

% Maksimalus mentes ilgis

w = (2*pi*RPM)/60; %[rad/s]

R_max = v_max / w; %[m]

% Sraigto disko plotas

W_cr = m_cr * g; %[N]

F = W_cr / (p * g); %[m^2]

% Mentes ilgis

R = sqrt(F/pi); %[m]

if R > R_max
    fprintf('Error R > R_max')
else
    fprintf('R < R_max == %4.2f m < %4.2f m \n', R, R_max)
end

% Oro srautas aplink mente

V = w * R; %[m/s]

fprintf('w = %4.2f m/s \n', w)
fprintf('R_max = %4.2f m \n', R_max)
fprintf('F = %4.2f m^2 \n', F)
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

P_gal_AG = P_gal / 0.7355; %[AG]

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
F_c = 0.5 * m_ment * w^2 * R^2;   % [N]

% Koordinate per mente nuo saknies iki galo
x = linspace(0, R, 500); %[m]

% Asine jega nuo iscentriniu jegu
Nx = 0.5 * m_ment * w^2 .* (R^2 - x.^2); %[N]

% Saknyje didziausia reiksme
Nx_root = Nx(1); %[N]

fprintf('Nx_root = %4.2f N \n', Nx_root)

figure;
plot(x, Nx, 'LineWidth', 2)
grid on
xlabel('Atstumas nuo mentes saknies x [m]')
ylabel('Asine jega N_x [N]')
title('Asines jegos pasiskirstymas menteje')

%% SKERSINIU JEGU IR LENKIMO MOMENTU SKAICIAVIMAS

rho = 1.225; %[kg/m^3]

% Pasirinktas keliamosios jegos koeficiento pasiskirstymas
Cl_root = 0.9;
Cl_tip  = 0.1;

% Tiesinis Cl kitimas per mente
Cl = Cl_root + (Cl_tip - Cl_root) * (x / R);

% Vietinis greitis kiekviename taske
Vx = w .* x; %[m/s]

% Zalias (nesukalibruotas) paskirstytas apkrovimas
q_raw = 0.5 * rho .* Vx.^2 .* C .* Cl; %[N/m]

% Kalibravimas
L_raw = trapz(x, q_raw);      %[N]
scale = T_ment / L_raw;

q = scale .* q_raw;           %[N/m]

% Patikra
L_blade = trapz(x, q);        %[N]

fprintf('L_raw = %4.2f N \n', L_raw)
fprintf('scale = %4.4f \n', scale)
fprintf('L_blade = %4.2f N \n', L_blade)

%% SKERSINE JEGA Q(x)
% Q(x) = integral nuo x iki R q(s) ds

Q = flip(cumtrapz(flip(x), flip(q))); %[N]
Q_root = Q(1);

fprintf('Q_root = %4.2f N \n', Q_root)

figure;
plot(x, Q, 'LineWidth', 2)
grid on
xlabel('Atstumas nuo mentes saknies x [m]')
ylabel('Skersine jega Q(x) [N]')
title('Skersines jegos pasiskirstymas menteje')

%% LENKIMO MOMENTAS M(x)
% M(x) = integral nuo x iki R Q(s) ds

M = flip(cumtrapz(flip(x), flip(Q))); %[N*m]
M_root = M(1);

fprintf('M_root = %4.2f N*m \n', M_root)

figure;
plot(x, M, 'LineWidth', 2)
grid on
xlabel('Atstumas nuo mentes saknies x [m]')
ylabel('Lenkimo momentas M(x) [N*m]')
title('Lenkimo momento pasiskirstymas menteje')
