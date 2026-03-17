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

% -----------------------------
% MANEVRINE GAUBIAMOJI
% -----------------------------

% Teigiama smukos kreive iki Va
V1 = linspace(0, Va, 300);
n1 = (V1 ./ Vs).^2;

% Virsutine horizontali atkarpa nuo Va iki Vd
V2 = linspace(Va, Vd, 300);
n2 = n_pos * ones(size(V2));

% Neigiama smukos kreive iki Vc_neg
V3 = linspace(0, Vc, 300);
n3 = n_neg * (V3 ./ Vc).^2;

% Neigiama tiesine atkarpa nuo Vc iki Vd, mazejanti iki 0
V4 = linspace(Vc, Vd, 300);
n4 = n_neg * (1 - (V4 - Vc) / (Vd - Vc));

% Pagrindiniai manevrines gaubiamosios taskai
A = [Vs, 1];
B = [Va, n_pos];
C = [Vc, n_pos];
D = [Vd, n_pos];
E = [Vc, n_neg];
F = [Vd, 0];

% -----------------------------
% GUST ENVELOPE
% -----------------------------
% Reikia:
% c_bar  - vidutine aerodinamine styga [m]
% a      - CLalpha [1/rad]
% Ude    - projektinis gusio greitis [m/s]

% Vidutine aerodinam. styga trapeciniam sparnui
c_bar = (2/3) * C_root * ((1 + (1/lam) + (1/lam)^2) / (1 + (1/lam)));

% CLalpha is XFLR5 arba artima prielaida
a = 3.22;              % [1/rad] <-- pasikeisk pagal savo XFLR5 jei turi

% Projektiniai gusio greiciai normal kategorijai
Ude_Vc = 15.24;       % [m/s] = 50 ft/s prie Vc
Ude_Vd = 7.62;        % [m/s] = 25 ft/s prie Vd

% Sparno apkrova SI sistema
WS_SI = W_N / S;      % [N/m^2]

% Mases santykio parametras mu_g
mu_g = (2 * WS_SI) / (rho * c_bar * a * g);

% Gusio lengvinimo koeficientas
K_g = (0.88 * mu_g) / (5.3 + mu_g);

% Gusio perkrovos prie Vc
delta_n_Vc = (K_g * rho * Vc * Ude_Vc * a) / (2 * WS_SI);
n_gust_pos_Vc = 1 + delta_n_Vc;
n_gust_neg_Vc = 1 - delta_n_Vc;

% Gusio perkrovos prie Vd
delta_n_Vd = (K_g * rho * Vd * Ude_Vd * a) / (2 * WS_SI);
n_gust_pos_Vd = 1 + delta_n_Vd;
n_gust_neg_Vd = 1 - delta_n_Vd;

% Tiesines gusiu kreives nuo 0 iki Vc ir nuo Vc iki Vd
Vg1 = linspace(0, Vc, 300);
n_gust_pos_1 = 1 + (n_gust_pos_Vc - 1) * (Vg1 / Vc);
n_gust_neg_1 = 1 + (n_gust_neg_Vc - 1) * (Vg1 / Vc);

Vg2 = linspace(Vc, Vd, 300);
n_gust_pos_2 = n_gust_pos_Vc + (n_gust_pos_Vd - n_gust_pos_Vc) * ((Vg2 - Vc) / (Vd - Vc));
n_gust_neg_2 = n_gust_neg_Vc + (n_gust_neg_Vd - n_gust_neg_Vc) * ((Vg2 - Vc) / (Vd - Vc));

% -----------------------------
% V_B skaiciavimas
% -----------------------------
% V_B = sankirta tarp teigiamos smukos kreives ir teigiamos gust kreives (0 iki Vc)

Vb_fun = @(V) (V./Vs).^2 - (1 + (n_gust_pos_Vc - 1) .* (V./Vc));

% ieskom sankirtos intervale [0, Vc]
Vb = fzero(Vb_fun, [0.1, Vc]);

% Perkrova ties V_B
n_B = (Vb / Vs)^2;

fprintf('Vb = %.3f m/s\n', Vb)
fprintf('n(Vb) = %.3f\n', n_B)

% -----------------------------
% GALUTINE ISORINE FLIGHT ENVELOPE RIBA
% -----------------------------
% Sukuriam bendra greicio asi
V_env = linspace(0, Vd, 1200);

% Manevrine virsutine riba
n_man_pos = zeros(size(V_env));
for i = 1:length(V_env)
    if V_env(i) <= Va
        n_man_pos(i) = (V_env(i)/Vs)^2;
    else
        n_man_pos(i) = n_pos;
    end
end

% Manevrine apatine riba
n_man_neg = zeros(size(V_env));
for i = 1:length(V_env)
    if V_env(i) <= Vc
        n_man_neg(i) = n_neg * (V_env(i)/Vc)^2;
    else
        n_man_neg(i) = n_neg * (1 - (V_env(i)-Vc)/(Vd-Vc));
    end
end

% Gust virsutine riba
n_gust_pos = zeros(size(V_env));
for i = 1:length(V_env)
    if V_env(i) <= Vc
        n_gust_pos(i) = 1 + (n_gust_pos_Vc - 1) * (V_env(i)/Vc);
    else
        n_gust_pos(i) = n_gust_pos_Vc + (n_gust_pos_Vd - n_gust_pos_Vc) * ((V_env(i)-Vc)/(Vd-Vc));
    end
end

% Gust apatine riba
n_gust_neg = zeros(size(V_env));
for i = 1:length(V_env)
    if V_env(i) <= Vc
        n_gust_neg(i) = 1 + (n_gust_neg_Vc - 1) * (V_env(i)/Vc);
    else
        n_gust_neg(i) = n_gust_neg_Vc + (n_gust_neg_Vd - n_gust_neg_Vc) * ((V_env(i)-Vc)/(Vd-Vc));
    end
end

% Galutine isorine riba:
% virsuje imam didesne, apacioje imam mazesne
n_env_pos = max(n_man_pos, n_gust_pos);
n_env_neg = min(n_man_neg, n_gust_neg);
% -----------------------------
% ISVEDIMAS
% -----------------------------
fprintf('\nGUST ENVELOPE:\n')
fprintf('----------------------------------------\n')
fprintf('c_bar = %.3f m\n', c_bar)
fprintf('a = %.3f 1/rad\n', a)
fprintf('mu_g = %.3f\n', mu_g)
fprintf('K_g = %.3f\n', K_g)
fprintf('n_gust_pos_Vc = %.3f\n', n_gust_pos_Vc)
fprintf('n_gust_neg_Vc = %.3f\n', n_gust_neg_Vc)
fprintf('n_gust_pos_Vd = %.3f\n', n_gust_pos_Vd)
fprintf('n_gust_neg_Vd = %.3f\n', n_gust_neg_Vd)

%% Braizymas 

figure;
hold on
grid on
box on

set(gca,'FontSize',12)
set(gca,'LineWidth',1.2)
set(gca,'GridLineStyle','-')
set(gca,'MinorGridLineStyle',':')
set(gca,'GridAlpha',0.25)
set(gca,'MinorGridAlpha',0.15)

% Simple 3-color scheme
maneuver_col = [0 0 1];   % blue
gust_col     = [1 0 0];   % red
env_col      = [0 0 0];   % black

% Manevrine gaubiamoji
h_man = plot(V1,n1,'Color',maneuver_col,'LineWidth',2);
plot(V2,n2,'Color',maneuver_col,'LineWidth',2)
plot(V3,n3,'Color',maneuver_col,'LineWidth',2)
plot(V4,n4,'Color',maneuver_col,'LineWidth',2)
plot([Vd Vd],[0 n_pos],'Color',maneuver_col,'LineWidth',2)

% Gust envelope
h_gust = plot(Vg1,n_gust_pos_1,'--','Color',gust_col,'LineWidth',2);
plot(Vg2,n_gust_pos_2,'--','Color',gust_col,'LineWidth',2)
plot(Vg1,n_gust_neg_1,'--','Color',gust_col,'LineWidth',2)
plot(Vg2,n_gust_neg_2,'--','Color',gust_col,'LineWidth',2)

% Galutine isorine flight envelope
h_env = plot(V_env,n_env_pos,'Color',env_col,'LineWidth',2.5);
plot(V_env,n_env_neg,'Color',env_col,'LineWidth',2.5)
plot([Vd Vd],[n_env_neg(end) n_env_pos(end)],'Color',env_col,'LineWidth',2.5)

% Taskai
plot(A(1),A(2),'ko','MarkerFaceColor','k')
plot(B(1),B(2),'ko','MarkerFaceColor','k')
plot(C(1),C(2),'ko','MarkerFaceColor','k')
plot(D(1),D(2),'ko','MarkerFaceColor','k')
plot(E(1),E(2),'ko','MarkerFaceColor','k')
plot(F(1),F(2),'ko','MarkerFaceColor','k')

% Vb taskas
plot(Vb,n_B,'ko','MarkerFaceColor','k','MarkerSize',8)

% Tasku pavadinimai
text(A(1), A(2), '  Vs', 'FontSize', 10)
text(B(1), B(2), '  Va', 'FontSize', 10)
text(C(1), C(2), '  Vc', 'FontSize', 10)
text(D(1), D(2), '  Vd', 'FontSize', 10)
text(E(1), E(2), '  Vc_{neg}', 'FontSize', 10)
text(F(1), F(2), '  0', 'FontSize', 10)
text(Vb, n_B, '  Vb', 'FontSize', 10)

% Pagalbines linijos
yline(0, 'k-')
yline(1, 'k--')
xline(Vs, 'k--')
xline(Va, 'k--')
xline(Vb, 'k--')
xline(Vc, 'k--')
xline(Vd, 'k--')

legend([h_man h_gust h_env], ...
       {'Manevrų gaubiamoji','Gūsinė gaubiamoji','Pilnoji skrydžio gaubiamoji'}, ...
       'Location','northwest')

xlabel('Greitis V [m/s]')
ylabel('Perkrovos koeficientas n')
title('Pilna perkrovos gaubiamoji')

xlim([0, 1.1*Vd])

y_min = min([n_env_neg, n_gust_neg, n_man_neg]) * 1.15;
y_max = max([n_env_pos, n_gust_pos, n_man_pos]) * 1.15;
ylim([y_min, y_max])

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

%% ------------------------------------------------------------
% SPYRIO ITAKA I LENKIMO MOMENTA - VIENAI SPARNO PUSEI
% ------------------------------------------------------------

idx_half = y >= 0;
y_half   = y(idx_half);
q_half   = q(idx_half);
BM_half  = BM(idx_half);

idx_out = y_half >= L_sp;

% spyrio vertikali reakcija - perimta isorines dalies apkrova
R_sp = trapz(y_half(idx_out), q_half(idx_out));

% spyrio asine jega
F_sp = R_sp / sin(alpha_rad);

% pakoreguotas momentas
BM_eff_half = BM_half;

idx_in = y_half <= L_sp;
BM_eff_half(idx_in) = BM_half(idx_in) - R_sp .* (L_sp - y_half(idx_in));

fprintf('\nSPYRIO ITAKA:\n');
fprintf('----------------------------------------\n');
fprintf('R_sp = %.3f N\n', R_sp);
fprintf('F_sp = %.3f N\n', F_sp);
fprintf('M_root be spyrio = %.3f Nm\n', BM_half(1));
fprintf('M_root su spyriu = %.3f Nm\n', BM_eff_half(1));

% grafikas
figure;
plot(y, q, 'LineWidth', 1.5);
grid on;
xlabel('Spanwise position y [m]');
ylabel('Distributed load q(y) [N/m]');
title('Wing load distribution');

% lenkimo momentas
figure;
plot(y_half, BM_half, 'LineWidth', 1.2); hold on;
plot(y_half, BM_eff_half, 'LineWidth', 1.8);
xline(0,'k--');
xline(L_sp,'k--');
grid on;
xlabel('Spanwise position y [m]');
ylabel('Bending moment M(y) [N m]');
title('Wing bending moment diagram');
legend('Be spyrio','Su spyriu','Location','best');

%% ============================================================
%LONZERONO SKAICIAVIMAS

H_spar = H_lonz;   % [m] sparo aukstis imamas is 1 kodo

% ------------------------------------------------------------
% 2) KRITINIS LENKIMO MOMENTAS IS XFLR5 DUOMENU
% ------------------------------------------------------------
% BM = data(:, col_BM);

% Imam didziausia absoliutine reiksme kaip kritine
[M_abs, idx_M] = max(abs(BM_eff_half));
M = M_abs;
y_crit = y_half(idx_M);


N = 0;              % [N] asine jega (siame variante ignoruojama)

fprintf('\nSPARO STIPRUMO SKAICIAVIMUI:\n');
fprintf('----------------------------------------\n');
fprintf('H_spar = %.6f m (is H_lonz)\n', H_spar);
fprintf('Kritinis lenkimo momentas M = %.3f Nm\n', M);
fprintf('Jis yra ties y = %.3f m\n', y_crit);

%% ------------------------------------------------------------
% 3) SKERSPJUvio PRIELAIDOS
% ------------------------------------------------------------

x1 = 0.20 * H_spar;      % lentynos plotis [m]
x2 = 0.02 * H_spar;      % sieneles storis [m]
x4 = 2 * x1;             % apsiuvos efektyvus plotis [m]

t_skin = 0.0015;         % apsiuvos storis [m]

% web aukstis + lentynos auksciai turi sudaryti visa sparo auksti
y2 = 0.60 * H_spar;      % sieneles aukstis [m]
y1 = (H_spar - y2) / 2;  % lentynos aukstis [m]

fprintf('\nSKERSPJUVIO MATMENYS:\n');
fprintf('----------------------------------------\n');
fprintf('Lentynos plotis x1 = %.6f m\n', x1);
fprintf('Sieneles storis x2 = %.6f m\n', x2);
fprintf('Apsiuvos plotis x4 = %.6f m\n', x4);
fprintf('Lentynos aukstis y1 = %.6f m\n', y1);
fprintf('Sieneles aukstis y2 = %.6f m\n', y2);
fprintf('Apsiuvos storis t_skin = %.6f m\n', t_skin);

% sluoksniai is apacios i virsu:
% 1 bottom skin, 2 bottom cap, 3 web, 4 top cap, 5 top skin
b = [x4,     x1,   x2,   x1,     x4];
t = [t_skin, y1,   y2,   y1, t_skin];
nLayers = numel(b);
A = b .* t;          % [m^2]

%% ------------------------------------------------------------
% 4) MEDZIAGOS
% ------------------------------------------------------------
E_7075 = 72e9;       % [Pa]
E_2024 = 73e9;       % [Pa]
E_6061 = 69e9;       % [Pa]

Re_7075 = 503e6;     % [Pa]
Re_2024 = 324e6;     % [Pa]
Re_6061 = 275e6;     % [Pa]

SF = 1.5;            % safety factor

%  bottom skin, bottom cap, web, top cap, top skin
E = [E_7075, E_2024, E_6061, E_2024, E_7075];

sigma_allow = [Re_7075, Re_2024, Re_6061, Re_2024, Re_7075] / SF;

%% ------------------------------------------------------------
% 5) BAZINIAI SKERSPJUVIO PARAMETRAI
% ------------------------------------------------------------
B_i = A .* E;        % [N]
B   = sum(B_i);      % [N]

H = sum(t);          % bendras aukstis
y_cent = zeros(1, nLayers);

y_running = 0;
for i = 1:nLayers
    y_cent(i) = y_running + t(i)/2;
    y_running = y_running + t(i);
end

%% ------------------------------------------------------------
% 6) NEUTRALIOJI ASIS
% ------------------------------------------------------------
y_NA = sum(B_i .* y_cent) / B;   % [m]

% atstumas nuo NA
y_iN = y_cent - y_NA;

%% ------------------------------------------------------------
% 7) LENKIMO STANDIS
% ------------------------------------------------------------
I_cent = (b .* t.^3) / 12;
I_i    = I_cent + A .* (y_iN.^2);

D_i = E .* I_i;
D   = sum(D_i);

%% ------------------------------------------------------------
% 8) ITEMPIAI
%    sigma_i = E_i * (N/B + M*y_i/D)
% ------------------------------------------------------------
epsilon_N = N / B;
sigma_i   = E .* (epsilon_N + (M .* y_iN) / D);   % [Pa]

% krastiniai pluostai
y_top = H - y_NA;
y_bot = -y_NA;

sigma_bot = E(1)   * (epsilon_N + M*y_bot/D);
sigma_top = E(end) * (epsilon_N + M*y_top/D);

%% ------------------------------------------------------------
% 9) ISNAUDOJIMAS
% ------------------------------------------------------------
util_sigma_layers = abs(sigma_i) ./ sigma_allow;
util_sigma_bot    = abs(sigma_bot) / sigma_allow(1);
util_sigma_top    = abs(sigma_top) / sigma_allow(end);

%% ------------------------------------------------------------
% 10) REZULTATU ISVEDIMAS
% ------------------------------------------------------------
fprintf('\n--- MULTILAYER SECTION STRENGTH CHECK (NORMAL STRESS) ---\n');
fprintf('Neutral axis from bottom y_NA = %.6f m\n', y_NA);
fprintf('Axial stiffness B             = %.3e N\n', B);
fprintf('Bending stiffness D           = %.3e Nm^2\n', D);
fprintf('Critical moment M             = %.3f Nm\n', M);
fprintf('Critical location y           = %.3f m\n\n', y_crit);

for i = 1:nLayers
    fprintf(['Layer %d: A = %.3e m^2, y_cent = %.6f m, ' ...
             'sigma = %+8.2f MPa, util = %.3f\n'], ...
        i, A(i), y_cent(i), sigma_i(i)/1e6, util_sigma_layers(i));
end

fprintf('\nExtreme fibres:\n');
fprintf('Bottom sigma = %+8.2f MPa, utilisation = %.3f\n', ...
    sigma_bot/1e6, util_sigma_bot);
fprintf('Top    sigma = %+8.2f MPa, utilisation = %.3f\n', ...
    sigma_top/1e6, util_sigma_top);

%% ------------------------------------------------------------
% 11) ITEMPIU PASISKIRSTYMO GRAFIKAS sigma(y)
% ------------------------------------------------------------
y_iface = [0, cumsum(t)];

nPts = 300;
y_vec = linspace(0, H, nPts);
sigma_y = zeros(size(y_vec));

for k = 1:nPts
    yy = y_vec(k);

    idx = find(yy >= y_iface(1:end-1) & yy <= y_iface(2:end), 1, 'first');

    % jei del apvalinimo paskutinis taskas nepapuola
    if isempty(idx)
        idx = nLayers;
    end

    y_rel = yy - y_NA;
    epsilon_y = epsilon_N + (M * y_rel) / D;
    sigma_y(k) = E(idx) * epsilon_y;
end

sigma_MPa = -sigma_y / 1e6;

figure;
hold on;
grid on;
plot(sigma_MPa, y_vec, 'LineWidth', 1.8);

xline(0, '--k', 'LineWidth', 1);
yline(y_NA, ':k', 'LineWidth', 1);

xlabel('\sigma [MPa]');
ylabel('y [m] (from bottom of spar)');
title('Normal stress distribution \sigma(y) in spar cross-section');

set(gca, 'YDir', 'normal');
hold off;
