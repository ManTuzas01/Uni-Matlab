%% Kintamieji parametrai
close all
clear
clc

% Duoti parametrai
m = 1700;           % [kg]
p_load = 70;        % max ikrovimas i sparna [kg/m^2]
n_e = 7;            % eksploatacinis perkrovos koeficientas
AR = 7;             % sparno proilgis
lam = 1.15;         % sparno susiaurejimas (C_root / C_tip)
H_sp = 1.0;         % spyrio aukstis [m]
L_sp = 1.5;         % spyrio ilgis [m]
g = 9.81;           % gravitacija [m/s^2]

t_prof = 11;        % sparno profilio storis [%]
f = 1.5;            % atsargos koeficientas

% Kuro bako duomenys
L_kuro = 1.500;     % [m]
H_kuro = 0.128067;  % [m]
w_kuro = 1.036687;  % [m]

%% 2.1 Sparno geometriniai parametrai (pagal pavyzdi)

% 1) Sparno plotas
S = m / p_load;                     % [m^2]

% 2) Sparno mostas
b = sqrt(S * AR);                   % [m]

% 3) Puses mostas
L = b / 2;                          % [m]

% 4) Galine styga
% Is trapecinio sparno ploto formules:
C_tip = (2 * S) / ((lam + 1) * b);  % [m]

% 5) Saknine styga
C_root = lam * C_tip;               % [m]

% 6) Styga ties spyriu

x = (C_root - C_tip) / 2;           % [m]
y = (L_sp * x) / L;                 % [m]
C_sp = C_tip + 2 * y;               % [m]

% Tas pats butu ir taip:
% C_sp = C_root - (C_root - C_tip)*(L_sp/L);

% 7) Maksimalus profilio storis ties spyriu
H_lonz = (t_prof / 100) * C_sp;     % [m]

% 8) Spyrio kampas
alpha_rad = atan(H_sp / L_sp);      % [rad]
alpha_deg = rad2deg(alpha_rad);     % [deg]



%% Kuro bako skaiciavimas

V_kuro = H_kuro * w_kuro * L_kuro;  % [m^3]
litrai = V_kuro * 1000;             % [L]

%% Rezultatu isvedimas

fprintf('\nAPSKAICIUOTA SPARNO GEOMETRIJA:\n')
fprintf('----------------------------------------\n')
fprintf('Sparno plotas S = %.3f m^2\n', S)
fprintf('Sparno mostas b = %.3f m\n', b)

fprintf('\nSTYGOS:\n')
fprintf('Saknine styga C_root = %.3f m\n', C_root)
fprintf('Styga ties spyriu C_sp = %.3f m\n', C_sp)
fprintf('Galine styga C_tip = %.3f m\n', C_tip)

fprintf('\nSPYRIS:\n')
fprintf('Spyrio kampas = %.3f rad\n', alpha_rad)
fprintf('Spyrio kampas = %.3f deg\n', alpha_deg)

fprintf('\nKURO BAKAS:\n')
fprintf('Kuro bako turis = %.3f m^3\n', V_kuro)
fprintf('Kuro kiekis = %.3f L\n', litrai)

%% Perkrovu ir greiciu skaiciavimas

rho = 1.225;                  % [kg/m^3]

% Trukstami dydziai
W_N = m * g;                  % [N]
W_kN = W_N / 1000;            % [kN]
n_sk = n_e * f;               % [-]
L_free = L - L_sp;            % [m]

% Vienetu konvertavimas i imperial (CS-23 formulems)
N_to_lbf = 0.224809;
m2_to_ft2 = 10.7639;
kt_to_ms = 0.514444;
ms_to_kt = 1.94384;

W_lbf = W_N * N_to_lbf;       % [lbf]
S_ft2 = S * m2_to_ft2;        % [ft^2]
WS_imp = W_lbf / S_ft2;       % [lbf/ft^2]

% Ribines perkrovos
n_pos = 2.1 + 24000 / (W_lbf + 10000);
n_pos = min(n_pos, 3.8);
n_neg = -0.4 * n_pos;

% Greiciai
Vc_kt = 33 * sqrt(WS_imp);    % [kt]
Vc = Vc_kt * kt_to_ms;        % [m/s]

Vd = 1.40 * Vc;               % [m/s]
Vd_kt = Vd * ms_to_kt;        % [kt]

% CLmax is XFLR5, apskaiciuok ir Vs bei Va
CLmax = 1.234;                 % <-- pakeisti i reiksme is XFLR5
Vs = sqrt((2 * W_N) / (rho * S * CLmax));   % [m/s]
Va = Vs * sqrt(n_pos);                       % [m/s]

% Isvedimas
fprintf('\nPERKROVOS IR GREICIAI:\n')
fprintf('----------------------------------------\n')
fprintf('W = %.3f kN\n', W_kN)
fprintf('W/S = %.3f lbf/ft^2\n', WS_imp)
fprintf('n_pos = %.3f\n', n_pos)
fprintf('n_neg = %.3f\n', n_neg)
fprintf('Vs = %.3f m/s\n', Vs)
fprintf('Va = %.3f m/s\n', Va)
fprintf('Vc = %.3f m/s (%.3f kt)\n', Vc, Vc_kt)
fprintf('Vd = %.3f m/s (%.3f kt)\n', Vd, Vd_kt)

%% V-n gaubiamoji

% Teigiama smukos kreive iki Va
V1 = linspace(0, Va, 300);
n1 = (V1 ./ Vs).^2;

% Virsutine horizontali atkarpa nuo Va iki Vd
V2 = linspace(Va, Vd, 300);
n2 = n_pos * ones(size(V2));

% Neigiama horizontali atkarpa iki Vc
V3 = linspace(0, Vc, 300);
n3 = n_neg * ones(size(V3));

% Neigiama tiesine atkarpa nuo Vc iki Vd, mazejanti iki 0
V4 = linspace(Vc, Vd, 300);
n4 = n_neg * (1 - (V4 - Vc) / (Vd - Vc));

% Pagrindiniai taskai
A = [Vs, 1];
B = [Va, n_pos];
C = [Vc, n_pos];
D = [Vd, n_pos];
E = [Vc, n_neg];
F = [Vd, 0];

% Flight envelope braizymas
figure;
hold on
grid on
box on

% Kreives
plot(V1, n1, 'b', 'LineWidth', 1.8)
plot(V2, n2, 'b', 'LineWidth', 1.8)
plot(V3, n3, 'b', 'LineWidth', 1.8)
plot(V4, n4, 'b', 'LineWidth', 1.8)

% Vertikali uzdarymo linija ties Vd
plot([Vd Vd], [0 n_pos], 'b', 'LineWidth', 1.8)

% Taskai
plot(A(1), A(2), 'ro', 'MarkerFaceColor', 'r')
plot(B(1), B(2), 'ro', 'MarkerFaceColor', 'r')
plot(C(1), C(2), 'ro', 'MarkerFaceColor', 'r')
plot(D(1), D(2), 'ro', 'MarkerFaceColor', 'r')
plot(E(1), E(2), 'ro', 'MarkerFaceColor', 'r')
plot(F(1), F(2), 'ro', 'MarkerFaceColor', 'r')

% Tasku pavadinimai
text(A(1), A(2), '  A', 'FontSize', 10)
text(B(1), B(2), '  B', 'FontSize', 10)
text(C(1), C(2), '  C', 'FontSize', 10)
text(D(1), D(2), '  D', 'FontSize', 10)
text(E(1), E(2), '  E', 'FontSize', 10)
text(F(1), F(2), '  F', 'FontSize', 10)

% Pagalbines linijos
yline(0, 'k-')
yline(1, '--k')
xline(Vs, '--k')
xline(Va, '--k')
xline(Vc, '--k')
xline(Vd, '--k')

xlabel('Greitis V [m/s]')
ylabel('Perkrovos koeficientas n')
title('Orlaivio manevrine V-n gaubiamoji')

xlim([0, 1.1*Vd])
ylim([1.2*n_neg, 1.2*n_pos])

hold off
%% Grafikai
% Plot distributed lift load q(y) over the wing from XFLR5 export
% Put your XFLR5 txt data into a file, e.g. 'xflr5_results.txt'

rho = 1.225;     % air density [kg/m^3]
V   = 123.7;     % flight speed [m/s]
qInf = 0.5 * rho * V^2;

filename = 'D:/MainWing_a=8.50_v=123.70ms.csv';

% perskaityti .csv failą
data = readmatrix(filename, 'FileType', 'text');

% stulpelių ID nr
col_y     = 1;   % y-span
col_chord = 2;   % Chord
col_Cl    = 4;   % Cl
col_BM    = 12;  % BM

y     = data(:, col_y);
chord = data(:, col_chord);
Cl    = data(:, col_Cl);
BM = data(:, col_BM);

% panaikinti visus kitus stulpelius
valid = ~isnan(y) & ~isnan(chord) & ~isnan(Cl) & ~isnan(BM);
y     = y(valid);
chord = chord(valid);
Cl    = Cl(valid);
BM    = BM(valid);

% išrikiuoti pagal y
[y, idx] = sort(y);
chord = chord(idx);
Cl    = Cl(idx);
BM    = BM(idx);

% q [N/m]
q = qInf .* chord .* Cl;

% grafikas
figure;
plot(y, q, 'LineWidth', 1.5);
grid on;
xlabel('Spanwise position y [m]');
ylabel('Distributed load q(y) [N/m]');
title('Wing load distribution');

% lenkimo momentas
figure;
plot(y, BM, 'LineWidth', 1.5);
grid on;
xlabel('Spanwise position y [m]');
ylabel('Bending moment M(y) [N m]');
title('Wing bending moment diagram');

% galutinis
figure;
subplot(2,1,1);
plot(y, q, 'LineWidth', 1.5);
grid on;
xlabel('Spanwise position y [m]');
ylabel('q(y) [N/m]');
title('Wing load distribution');

subplot(2,1,2);
plot(y, BM, 'LineWidth', 1.5);
grid on;
xlabel('Spanwise position y [m]');
ylabel('M(y) [N m]');
title('Wing bending moment diagram');
