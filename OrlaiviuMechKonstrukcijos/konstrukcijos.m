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

V_kuro = 0.1;  % [m^3]
litrai = V_kuro * 1000;             % [L]

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

CL_max = 1.36825;

Vs = sqrt((2 * W_N) / (rho * S * CL_max));

n_pos = n_e;
n_neg = -n_e;
Va = Vs * sqrt(n_pos); 

Vc_kt = 33 * sqrt(WS_imp);     % minimum design cruise speed in knots
Vc = Vc_kt * kt_to_ms;         % [m/s]

Vd = 1.40 * Vc;                % [m/s]
Vd_kt = Vd * ms_to_kt;         % [knots]

L_Vd = n_pos * W_N;                          % [N]
CL_Vd = (2 * L_Vd) / (rho * Vd^2 * S);      % [-]

alpha_cl = 12.5;

%% Skrydžio perkrovu gaubiamoji

% Simetrine riba
n_neg = -n_pos;

% Smukos kreive iki Va
V_stall = linspace(0, Va, 300);
n_stall_pos = (V_stall.^2 ./ Vs^2);
n_stall_neg = -(V_stall.^2 ./ Vs^2);

% Horizontalios ribos nuo Va iki Vd
V_lim = linspace(Va, Vd, 300);
n_lim_pos = n_pos * ones(size(V_lim));
n_lim_neg = n_neg * ones(size(V_lim));

% Pagrindiniai taskai
A = [Vs, 1];
B = [Va, n_pos];
C = [Vc, n_pos];
D = [Vd, n_pos];

A2 = [Vs, -1];
B2 = [Va, n_neg];
C2 = [Vc, n_neg];
D2 = [Vd, n_neg];

%% XFLR5 duomenu nuskaitymas ir braizymas

filename_xflr = 'D:/MainWing_a=12.50_v=90.00ms - Copy.txt';

data = readmatrix(filename_xflr);

% Stulpeliai pagal tavo sutvarkyta txt faila
y_x     = data(:,1);    % y-span [m]
chord_x = data(:,2);    % chord [m]
Cl_x    = data(:,4);    % local Cl [-]

% Dinaminis slegis
V_xflr = 90.0;                  % [m/s]
q_inf = 0.5 * rho * V_xflr^2;   % [Pa]

% Paskirstyta apkrova
q_span = q_inf .* chord_x .* Cl_x;   % [N/m]


%% Kirpimo jega ir lenkimo momentas is q(y)

% Isrikiuojam pagal sparno koordinate
[y_sort, idx] = sort(y_x);
q_sort = q_span(idx);

% Imam tik desine puse nuo saknies iki galo
mask_right = y_sort >= 0;
y_right_raw = y_sort(mask_right);
q_right_raw = q_sort(mask_right);

% Panaikinam pasikartojancias y reiksmes, jei tokiu yra
[y_right_raw, ia] = unique(y_right_raw, 'stable');
q_right_raw = q_right_raw(ia);

if y_right_raw(1) > 0
    y_right_raw = [0; y_right_raw];
    q_right_raw = [q_right_raw(1); q_right_raw];
end

% Sudarom tankesni tinkleli glotnesnei kreivei
y_right = linspace(min(y_right_raw), max(y_right_raw), 400);

% Interpoliuojam apkrova
q_right = interp1(y_right_raw, q_right_raw, y_right, 'pchip');

% Kirpimo jega:
% ties sparno galu V = 0, todel integruojam nuo galo link saknies
V_right = flip(cumtrapz(flip(y_right), flip(q_right)));

% Lenkimo momentas:
% ties sparno galu M = 0, todel vel integruojam nuo galo link saknies
M_right = flip(cumtrapz(flip(y_right), flip(V_right)));


%% Spyrio apkrova ir pakoreguotas lenkimo momentas

% Pradinis lenkimo momentas ties saknimi be spyrio
M_root_0 = M_right(1);                       % [N·m]

% Reikalinga spyrio vertikali reakcija, kad momentas ties fiuzeliazu butu 0
F_sp_v = M_root_0 / L_sp;                    % [N]

% Spyrio asine jega
N_sp = F_sp_v / sin(alpha_rad);              % [N]

% Spyrio horizontali dedamoji
F_sp_h = N_sp * cos(alpha_rad);              % [N]

% Pakoreguotos kirpimo jegos ir momentai
V_right_sp = V_right;
M_right_sp = M_right;

% Ties taskais nuo saknies iki spyrio tvirtinimo tasko
idx_sp = y_right <= L_sp;

% Spyris mazina kirpimo jega ir lenkimo momenta vidineje dalyje
V_right_sp(idx_sp) = V_right_sp(idx_sp) - F_sp_v;
M_right_sp(idx_sp) = M_right_sp(idx_sp) - F_sp_v .* (L_sp - y_right(idx_sp));

%% Rezultatu isvedimas

%% pirma dalis
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

%% Antra dalis
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

%% Trecia dalis
fprintf('\nSPYRIO APKROVA:\n')
fprintf('----------------------------------------\n')
fprintf('Pradinis momentas ties saknimi be spyrio = %.3f N·m\n', M_root_0)
fprintf('Spyrio vertikali dedamoji F_sp_v = %.3f N\n', F_sp_v)
fprintf('Spyrio asine jega N_sp = %.3f N\n', N_sp)
fprintf('Spyrio horizontali dedamoji F_sp_h = %.3f N\n', F_sp_h)
fprintf('Pakoreguotas momentas ties saknimi = %.3f N·m\n', M_right_sp(1))

%% Skrydzio gaubiamosios braizymas

figure;
hold on
grid on
box on

% Virsutine dalis
plot(V_stall, n_stall_pos, 'b', 'LineWidth', 2)
plot(V_lim,   n_lim_pos,   'b', 'LineWidth', 2)

% Apatine dalis (simetriska)
plot(V_stall, n_stall_neg, 'b', 'LineWidth', 2)
plot(V_lim,   n_lim_neg,   'b', 'LineWidth', 2)

% Uzsidarymas ties Vd
plot([Vd Vd], [n_neg n_pos], 'b', 'LineWidth', 2)

% Pagalbines vertikalios linijos
xline(Vs, '--k')
xline(Va, '--k')
xline(Vc, '--k')
xline(Vd, '--k')

% X asis
yline(0, 'k-')

% Taskai
plot(A(1),  A(2),  'ro', 'MarkerFaceColor', 'r')
plot(B(1),  B(2),  'ro', 'MarkerFaceColor', 'r')
plot(C(1),  C(2),  'ro', 'MarkerFaceColor', 'r')
plot(D(1),  D(2),  'ro', 'MarkerFaceColor', 'r')

plot(A2(1), A2(2), 'ro', 'MarkerFaceColor', 'r')
plot(B2(1), B2(2), 'ro', 'MarkerFaceColor', 'r')
plot(C2(1), C2(2), 'ro', 'MarkerFaceColor', 'r')
plot(D2(1), D2(2), 'ro', 'MarkerFaceColor', 'r')

% Uzrasai
text(A(1),  A(2),  '  Vs')
text(B(1),  B(2),  '  Va')
text(C(1),  C(2),  '  Vc')
text(D(1),  D(2),  '  Vd')

xlabel('Greitis V [m/s]')
ylabel('Perkrovos koeficientas n')
title('Skrydžio perkrovų gaubiamoji')

%% apkrovos

figure;
plot(y_x, q_span, 'o-','LineWidth',1.2)
grid on
xlabel('Sparno koordinatė y [m]')
ylabel('Paskirstyta apkrova q(y) [N/m]')
title('Sparno apkrovos pasiskirstymas')

%% Kirpimo jegos ir lenkimo momento braizymas

figure;
plot(y_right, V_right, 'LineWidth', 1.5)
grid on
xlabel('Pusės sparno koordinatė y [m]')
ylabel('Kirpimo jėga V(y) [N]')
title('Sparno kirpimo jėgos diagrama')

figure;
plot(y_right, M_right, 'LineWidth', 1.5)
grid on
xlabel('Pusės sparno koordinatė y [m]')
ylabel('Lenkimo momentas M(y) [N·m]')
title('Sparno lenkimo momento diagrama')

%% Pakoreguotos diagramos su spyriu

figure;
plot(y_right, V_right, 'LineWidth', 1.2)
hold on
plot(y_right, V_right_sp, 'LineWidth', 1.5)
xline(L_sp, '--k')
grid on
xlabel('Pusės sparno koordinatė y [m]')
ylabel('Kirpimo jėga V(y) [N]')
title('Kirpimo jėga prieš ir po spyrio')
legend('Be spyrio','Su spyriu','Spyrio vieta')

figure;
plot(y_right, M_right, 'LineWidth', 1.2)
hold on
plot(y_right, M_right_sp, 'LineWidth', 1.5)
xline(L_sp, '--k')
grid on
xlabel('Pusės sparno koordinatė y [m]')
ylabel('Lenkimo momentas M(y) [N·m]')
title('Lenkimo momentas prieš ir po spyrio')
legend('Be spyrio','Su spyriu','Spyrio vieta')
