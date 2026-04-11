R < R_max == 5.35 m < 6.29 m 
w = 41.36 m/s 
R_max = 6.29 m 
F = 89.93 m^2 
R = 5.35 m 
V = 221.31 m/s 
T = 24525.00 N 
T_ment = 8175.00 N 
v_i = 10.55 m/s 
P_ideal = 258.75 kW 
P_prof = 64.69 kW 
P = 323.44 kW 
P_vel = 359.38 kW 
P_gal = 431.25 kW 
P_gal_AG = 586.34 AG 
Nx_root = 269376.35 N 
L_raw = 5296.50 N 
scale = 1.5435 
L_blade = 8175.00 N 
Q_root = 8175.00 N 
M_root = 28429.74 N*m 

--- SUJUNGIMO APKROVOS ---
Nx_root = 269376.35 N
Q_root  =  8175.00 N
Vienam varztui tenkanti rezultuojanti apkrova F_bolt = 134750.18 N

--- VARZTO PARINKIMAS ---
Minimalus varzto skersmuo is kirpimo =  16.91 mm
Parinktas standartinis varzto skersmuo =  18.00 mm
Skyles skersmuo =  19.00 mm

--- VARZTO KIRPIMO PATIKRA ---
Tikroji varzto kirpimo itempa =   264.77 MPa
Leistina varzto kirpimo itempa =   300.00 MPa

--- BEARING PATIKRA / REIKALINGI STORIAI ---
Minimalus centrines saknies dalies storis t_root_req =  28.37 mm
Pasirinktas centrines dalies storis t_root =  20.00 mm
Minimalus vienos auseles storis t_lug_req =  28.37 mm
Pasirinktas vienos auseles storis t_lug =  15.00 mm
Bearing itempa centrineje dalyje =   354.61 MPa
Bearing itempa auselese          =   472.81 MPa
Leistina bearing itempa         =   250.00 MPa

--- NETTO PJUVIO TEMPIMO PATIKRA ---
Minimalus saknines dalies plotis b_root_req =  52.69 mm

--- KRASHTO ISPLESIMO (SHEAR-OUT) PATIKRA ---
Minimalus krasto atstumas centrinei daliai e_root_req =  22.46 mm
Minimalus krasto atstumas auselese e_lug_req          =  29.94 mm

---  JUNGTIES MATMENYS ---
 varzto skersmuo d        =  18.00 mm
 skyles skersmuo d0       =  19.00 mm
 krasto atstumas e        =  36.00 mm
 atstumas tarp skyliu s   =  63.00 mm
 saknines dalies plotis b =  72.00 mm
 auseles aukstis h        =  91.00 mm
 jungties ilgis L         = 135.00 mm
Pasirinktas centrines dalies storis      =  20.00 mm
Pasirinktas vienos auseles storis        =  15.00 mm

--- SAUGOS KOEFICIENTAI ---
SF varzto kirpimui              =   1.13
SF bearing centrineje dalyje   =   0.71
SF bearing auselese            =   0.53
SF netto tempimui              =   1.57
SF shear-out centrineje dalyje =   1.60
SF shear-out auselese          =   1.20
>> clc

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
fprintf('V = %4.2f m/s \n', V)

%% KABOJIMUI REIKALINGA GALIA IR JEGA

rho = 1.225; %[kg/m^3]

% Traukos jega
T = W_cr; %[N]
T_ment = T / n_ment; %[N]

% Indukuotas greitis kabojime
v_i = sqrt(T / (2 * rho * F)); %[m/s]

% Ideali galia
P_ideal = T * v_i / 1000; %[kW]

% Profilio galia = 25% idealios galios
P_prof = 0.25 * P_ideal; %[kW]

% Bendra rotoriui reikalinga galia
P = P_ideal + P_prof; %[kW]

% Ivertinus transmisijos NK
P_vel = P / n_coef; %[kW]

% Su 20% rezervu
P_gal = 1.2 * P_vel; %[kW]

P_gal_AG = P_gal / 0.7355; %[AG]

fprintf('T = %4.2f N \n', T)
fprintf('T_ment = %4.2f N \n', T_ment)
fprintf('v_i = %4.2f m/s \n', v_i)
fprintf('P_ideal = %4.2f kW \n', P_ideal)
fprintf('P_prof = %4.2f kW \n', P_prof)
fprintf('P = %4.2f kW \n', P)
fprintf('P_vel = %4.2f kW \n', P_vel)
fprintf('P_gal = %4.2f kW \n', P_gal)
fprintf('P_gal_AG = %4.2f AG \n', P_gal_AG)
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

Q = -flip(cumtrapz(flip(x), flip(q))); %[N]
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

M = -flip(cumtrapz(flip(x), flip(Q))); %[N*m]
M_root = M(1);

fprintf('M_root = %4.2f N*m \n', M_root)

figure;
plot(x, M, 'LineWidth', 2)
grid on
xlabel('Atstumas nuo mentes saknies x [m]')
ylabel('Lenkimo momentas M(x) [N*m]')
title('Lenkimo momento pasiskirstymas menteje')
