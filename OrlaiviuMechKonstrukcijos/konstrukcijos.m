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

%% V-n gaubiamoji + gust envelope

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

% -----------------------------
% BRAIZYMAS
% -----------------------------
%% DARK MODE GRAFIKAS

bg = [0.08 0.08 0.08];     % tamsus fonas
axis_col = [0.9 0.9 0.9];  % ašių ir teksto spalva
grid_col = [0.35 0.35 0.35];

set(groot,'defaultFigureColor',bg)
set(groot,'defaultAxesColor',bg)

set(groot,'defaultAxesXColor',axis_col)
set(groot,'defaultAxesYColor',axis_col)

set(groot,'defaultTextColor',axis_col)

set(groot,'defaultAxesGridColor',grid_col)
set(groot,'defaultAxesMinorGridColor',grid_col)

set(groot,'defaultAxesFontSize',12)
set(groot,'defaultAxesLineWidth',1.2)

set(groot,'defaultLegendTextColor',axis_col)
set(groot,'defaultLegendColor',bg)
set(groot,'defaultLegendEdgeColor',axis_col)

figure;
hold on
grid on
box on
set(gca,'GridAlpha',0.35)
set(gca,'MinorGridAlpha',0.25)
set(gca,'GridLineStyle','-')
set(gca,'MinorGridLineStyle',':')

maneuver_col = [0 0.85 1];     % cyan
gust_col     = [1 0.55 0];     % orange
env_col      = [1 1 1];        % white
point_col    = [1 1 0];        % yellow

% Manevrine gaubiamoji
h_man = plot(V1,n1,'Color',maneuver_col,'LineWidth',2);
plot(V2,n2,'Color',maneuver_col,'LineWidth',2)
plot(V3,n3,'Color',maneuver_col,'LineWidth',2)
plot(V4,n4,'Color',maneuver_col,'LineWidth',2)
plot([Vd Vd],[0 n_pos],'Color',maneuver_col,'LineWidth',2)

% Gust envelope - raudona bruksniuota
h_gust = plot(Vg1,n_gust_pos_1,'--','Color',gust_col,'LineWidth',2);
plot(Vg2,n_gust_pos_2,'--','Color',gust_col,'LineWidth',2)
plot(Vg1,n_gust_neg_1,'--','Color',gust_col,'LineWidth',2)
plot(Vg2,n_gust_neg_2,'--','Color',gust_col,'LineWidth',2)

% Galutine isorine flight envelope - juoda stora
h_env = plot(V_env,n_env_pos,'Color',env_col,'LineWidth',3);
plot(V_env,n_env_neg,'Color',env_col,'LineWidth',3)
plot([Vd Vd],[n_env_neg(end) n_env_pos(end)],'Color',env_col,'LineWidth',3)

% Taskai
plot(A(1),A(2),'o','Color',point_col,'MarkerFaceColor',point_col)
plot(B(1),B(2),'o','Color',point_col,'MarkerFaceColor',point_col)
plot(C(1),C(2),'o','Color',point_col,'MarkerFaceColor',point_col)
plot(D(1),D(2),'o','Color',point_col,'MarkerFaceColor',point_col)
plot(E(1),E(2),'o','Color',point_col,'MarkerFaceColor',point_col)
plot(F(1),F(2),'o','Color',point_col,'MarkerFaceColor',point_col)

% Vb taskas
plot(Vb,n_B,'o','Color',[1 0 1],'MarkerFaceColor',[1 0 1],'MarkerSize',8)

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
yline(1, '--k')
xline(Vs, '--k')
xline(Va, '--k')
xline(Vb, '--m')
xline(Vc, '--k')
xline(Vd, '--k')


legend([h_man h_gust h_env], ...
       {'Manevrų gaubiamoji','Gūsinė gaubiamoji','Pilnoji skrydio gaubiamoji'}, ...
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
