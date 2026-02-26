%% Kintamieji parametrai
close all
clear
clc

% Duoti parametrai
m = 1700;           % [kg]
p_load = 70;        % sparno apkrova [kg/m^2]
n_e = 7;            % eksploatacinis perkrovos koeficientas
AR = 7;             % sparno pailgėjimas
lam = 1.15;         % siaurėjimas (C_root / C_tip)
H_sp = 1;           % spyrio aukštis [m]
L_sp = 1.5;         % spyrio ilgis [m]
g = 9.81;           % gravitacija
V_fuel = 200e-3;    % kuro tūris [m^3]

%% 2.1 Sparno geometriniai parametrai

% Svoris
W = m * g;                          % [N]

% Sparno plotas
p_load_N = p_load * g;              % [N/m^2]
S = W / p_load_N;                   % [m^2]

% Sparno mostas
b = sqrt(S * AR);                   % [m]

% Pusės mostas
L = b / 2;                          % [m]

% Galinė styga
C_tip = (2*S) / (b*(lam+1));       % [m]

% Šakninė styga
C_root = lam * C_tip;              % [m]

% Spyrio padėtis
b1 = L_sp;                          % [m]
b2 = L - b1;                        % laisvas galas [m]

% Styga ties spyriu (linijinė interpolacija)
C_sp = C_root - (C_root - C_tip)*(b1/L);

% Vidutinė aerodinaminė styga (MAC)
C_MAC = (2/3)*C_root*((1+lam+lam^2)/(1+lam));

%% Spyrio kampas

theta_sp = atan(H_sp / L_sp);      % [rad]
theta_sp_deg = rad2deg(theta_sp);  % [deg]

%% Sparno aukštis ir kuro bakas

% spar aukštis (14% nuo MAC)
H_spar = 0.14 * C_MAC;


%NENAUDOTI ŠITO, MODELIUOTI KURO BAKĄ PER SOLIDWORKS RANKA NE SKAIČIAIS
% kuro bako aukštis
H_tank = H_spar - 0.0108;

% kuro bako styga
C_tank = V_fuel / (b1 * H_tank);

%% Rezultatų išvedimas

fprintf('\nSPARNO PARAMETRAI:\n')

fprintf('Plotas S = %.3f m^2\n', S)
fprintf('Mostas b = %.3f m\n', b)

fprintf('\nStygos:\n')
fprintf('Šakninė C_root = %.3f m\n', C_root)
fprintf('Galinė C_tip = %.3f m\n', C_tip)
fprintf('Styga ties spyriu C_sp = %.3f m\n', C_sp)
fprintf('MAC = %.3f m\n', C_MAC)

fprintf('\nGeometrija:\n')
fprintf('Spyrio kampas = %.3f deg\n', theta_sp_deg)
fprintf('Laisvas galas b2 = %.3f m\n', b2)

fprintf('\nKuro bakas:\n')
fprintf('Aukštis = %.3f m\n', H_tank)
fprintf('Styga = %.3f m\n', C_tank)