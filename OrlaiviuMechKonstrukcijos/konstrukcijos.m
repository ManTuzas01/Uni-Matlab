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

V_kuro = L_kuro * H_kuro * w_kuro;  % [m^3]
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
A_pt = [Vs, 1];
B_pt = [Va, n_pos];
C_pt = [Vc, n_pos];
D_pt = [Vd, n_pos];

A2_pt = [Vs, -1];
B2_pt = [Va, n_neg];
C2_pt = [Vc, n_neg];
D2_pt = [Vd, n_neg];

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


%% =========================================================
%  ITEMPIŲ SKAIČIAVIMAS KRITINIAME PJŪVYJE (TIES SPYRIU)
% =========================================================

% Kritinis taškas – spyrio tvirtinimo vieta
y_crit = L_sp;   % [m]

% Vidinės jėgos kritiniame pjūvyje
V_crit = interp1(y_right, V_right_sp, y_crit, 'linear');
M_crit = interp1(y_right, M_right_sp, y_crit, 'linear');

% Jei nori ignoruoti ašinę jėgą sparo pjūvyje:
N_crit = 0;

% Jei nori įtraukti spyrio horizontalią dedamąją kaip ašinę jėgą sparui,
% gali vietoje viršutinės eilutės naudoti, pvz.:
% N_crit = -F_sp_h;
%
% Minusas reiškia gniuždymą, jei laikome tempimą teigiamu.
% Kol neaiški tiksli jėgų schema, saugiausia pradžiai palikti N_crit = 0.

%% ------------------------------------------------------------------------
%  GEOMETRY OF MULTILAYER SPAR SECTION
% -------------------------------------------------------------------------
H_spar = 0.21225;     % total spar height [m]

x1 = 0.2;             % spar cap width [m]
x2 = 0.03;            % web thickness [m]
x4 = 2*x1;            % skin width [m]

t_skin = 0.0015;      % skin thickness (top & bottom) [m]
y2_sec = 0.13;        % web height [m]
y1_sec = (H_spar - y2_sec)/2;   % spar cap thickness [m]

fprintf('\nKRITINIO PJŪVIO GEOMETRIJA:\n')
fprintf('Lentynos plotis x1 = %.6f m\n', x1)
fprintf('Sienelės plotis x2 = %.6f m\n', x2)
fprintf('Apsiuvos plotis x4 = %.6f m\n', x4)
fprintf('Lentynos aukštis y1 = %.6f m\n', y1_sec)
fprintf('Sienelės aukštis y2 = %.6f m\n', y2_sec)
fprintf('Apsiuvos storis t_skin = %.6f m\n', t_skin)

% sluoksniai iš apačios į viršų:
% 1 apatinė apsiūva, 2 apatinė lentyna, 3 sienelė, 4 viršutinė lentyna, 5 viršutinė apsiūva
b_sec = [x4,     x1,      x2,      x1,     x4];
t_sec = [t_skin, y1_sec,  y2_sec,  y1_sec, t_skin];
nLayers = numel(b_sec);

A_sec = b_sec .* t_sec;   % [m^2]

%% ------------------------------------------------------------------------
%  MATERIAL DATA
% -------------------------------------------------------------------------
E_7075 = 72e9;       % [Pa]
E_2024 = 73e9;       % [Pa]
E_6061 = 69e9;       % [Pa]

Re_7075 = 503e6;     % [Pa]
Re_2024 = 324e6;     % [Pa]
Re_6061 = 275e6;     % [Pa]

SF = 1.5;

E_sec = [E_7075, E_2024, E_6061, E_2024, E_7075];
sigma_allow = [Re_7075, Re_2024, Re_6061, Re_2024, Re_7075] / SF;

%% ------------------------------------------------------------------------
%  BASIC SECTION PROPERTIES
% -------------------------------------------------------------------------
B_i = A_sec .* E_sec;
B_ax = sum(B_i);

H_sec = sum(t_sec);
y_cent = zeros(1, nLayers);

y_running = 0;
for i = 1:nLayers
    y_cent(i) = y_running + t_sec(i)/2;
    y_running = y_running + t_sec(i);
end

%% ------------------------------------------------------------------------
%  NEUTRAL AXIS
% -------------------------------------------------------------------------
y_NA = sum(B_i .* y_cent) / B_ax;

% atstumai nuo neutraliosios ašies
y_iN = y_cent - y_NA;

%% ------------------------------------------------------------------------
%  SECOND MOMENT OF AREA AND BENDING STIFFNESS
% -------------------------------------------------------------------------
I_cent = (b_sec .* t_sec.^3) / 12;
I_i    = I_cent + A_sec .* (y_iN.^2);

D_i = E_sec .* I_i;
D   = sum(D_i);

%% ------------------------------------------------------------------------
%  NORMAL STRESSES IN CRITICAL SECTION
%  sigma_i = E_i * ( N/B + M*y_i/D )
% -------------------------------------------------------------------------
epsilon_N = N_crit / B_ax;
sigma_i   = E_sec .* (epsilon_N + (M_crit .* y_iN) / D);

% kraštiniai taškai
y_top = H_sec - y_NA;
y_bot = -y_NA;

sigma_bot = E_sec(1)   * (epsilon_N + M_crit * y_bot / D);
sigma_top = E_sec(end) * (epsilon_N + M_crit * y_top / D);

%% ------------------------------------------------------------------------
%  UTILISATION
% -------------------------------------------------------------------------
util_sigma_layers = abs(sigma_i) ./ sigma_allow;
util_sigma_bot    = abs(sigma_bot) / sigma_allow(1);
util_sigma_top    = abs(sigma_top) / sigma_allow(end);

%% ------------------------------------------------------------------------
%  RESULTS
% -------------------------------------------------------------------------
fprintf('\nKRITINIS PJŪVIS TIES SPYRIU:\n')
fprintf('----------------------------------------\n')
fprintf('Pjūvio vieta y = %.3f m\n', y_crit)
fprintf('Kirpimo jėga V_crit = %.3f N\n', V_crit)
fprintf('Lenkimo momentas M_crit = %.3f N·m\n', M_crit)
fprintf('Ašinė jėga N_crit = %.3f N\n', N_crit)

fprintf('\nDAUGIASLUOKSNIO PJŪVIO STIPRUMO PATIKRA:\n')
fprintf('Neutralioji ašis nuo apačios y_NA = %.6f m\n', y_NA)
fprintf('Ašinis standumas B = %.3e N\n', B_ax)
fprintf('Lenkimo standumas D = %.3e N·m^2\n\n', D)

for i = 1:nLayers
    fprintf(['Sluoksnis %d: A = %.3e m^2, y_cent = %.4f m, ' ...
             'sigma = %+8.2f MPa, util = %.3f\n'], ...
        i, A_sec(i), y_cent(i), sigma_i(i)/1e6, util_sigma_layers(i));
end

fprintf('\nKraštiniai pluoštai:\n');
fprintf('Apatinis sigma = %+8.2f MPa, utilisation = %.3f\n', ...
    sigma_bot/1e6, util_sigma_bot);
fprintf('Viršutinis sigma = %+8.2f MPa, utilisation = %.3f\n', ...
    sigma_top/1e6, util_sigma_top);

% didžiausia absoliutinė įtampa
sigma_max = max(abs([sigma_i, sigma_bot, sigma_top]));
util_max  = max([util_sigma_layers, util_sigma_bot, util_sigma_top]);

fprintf('\nDidžiausia |sigma| = %.2f MPa\n', sigma_max/1e6)
fprintf('Didžiausias utilisation = %.3f\n', util_max)

%% ------------------------------------------------------------------------
%  NORMAL STRESS DISTRIBUTION PLOT sigma(y)
% -------------------------------------------------------------------------
y_iface = [0, cumsum(t_sec)];

nPts = 300;
y_vec = linspace(0, H_sec, nPts);
sigma_y = zeros(size(y_vec));

for k = 1:nPts
    yloc = y_vec(k);

    idx = find(yloc >= y_iface(1:end-1) & yloc <= y_iface(2:end), 1, 'first');
    y_rel = yloc - y_NA;

    epsilon_y = epsilon_N + (M_crit * y_rel) / D;
    sigma_y(k) = E_sec(idx) * epsilon_y;
end

sigma_MPa = -sigma_y / 1e6;   % pagal tavo pasirinktą ženklų konvenciją

figure;
hold on
grid on
plot(sigma_MPa, y_vec, 'LineWidth', 1.8)

xline(0, '--k', 'LineWidth', 1)
yline(y_NA, ':k', 'LineWidth', 1)

xlabel('\sigma [MPa]')
ylabel('y [m] (nuo sparo apačios)')
title('Normalinių įtempių pasiskirstymas kritiniame pjūvyje ties spyriu')

set(gca, 'YDir', 'normal')

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
plot(A_pt(1),  A_pt(2),  'ro', 'MarkerFaceColor', 'r')
plot(B_pt(1),  B_pt(2),  'ro', 'MarkerFaceColor', 'r')
plot(C_pt(1),  C_pt(2),  'ro', 'MarkerFaceColor', 'r')
plot(D_pt(1),  D_pt(2),  'ro', 'MarkerFaceColor', 'r')

plot(A2_pt(1), A2_pt(2), 'ro', 'MarkerFaceColor', 'r')
plot(B2_pt(1), B2_pt(2), 'ro', 'MarkerFaceColor', 'r')
plot(C2_pt(1), C2_pt(2), 'ro', 'MarkerFaceColor', 'r')
plot(D2_pt(1), D2_pt(2), 'ro', 'MarkerFaceColor', 'r')

text(A_pt(1),  A_pt(2),  '  Vs')
text(B_pt(1),  B_pt(2),  '  Va')
text(C_pt(1),  C_pt(2),  '  Vc')
text(D_pt(1),  D_pt(2),  '  Vd')

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
