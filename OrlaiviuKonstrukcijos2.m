clear all;
close all;
clc;

%% ------------------------------------------------------------------------
%  GEOMETRY
% -------------------------------------------------------------------------
H_spar = 0.1845;     % total spar height [m]

x1 = 0.2;           % spar cap width [m]
x2 = 0.03;           % web thickness [m]
x4 = 2*x1;           % skin width [m]

t_skin = 0.0015;      % skin thickness (top & bottom) [m]
y2     = 0.13;  % web height [m]
y1     = (H_spar - y2)/2;       % spar cap thickness [m]

fprintf('Lentynos plotis x1 = %.6f m', x1)
fprintf('Sieneles plotis x2 = %.6f m', x2)
fprintf('Apsiuvos plotis x3 = %.6f m', x4)
fprintf('Lentynos aukstis y1 = %.6f m', y1)
fprintf('Sieneles aukstis y2 = %.6f m', y2)
fprintf('Apsiuvos aukstis y3 = %.6f m', t_skin)

% layer widths (b_i) and thicknesses (t_i) from BOTTOM to TOP:
% 1 bottom skin, 2 bottom cap, 3 web, 4 top cap, 5 top skin
b = [x4,     x1,   x2,   x1,     x4];
t = [t_skin, y1,   y2,   y1, t_skin];
nLayers = numel(b);
A = b .* t;          % [m^2]

%% ------------------------------------------------------------------------
%  MATERIAL DATA (PER LAYER)
% -------------------------------------------------------------------------
E_7075 = 72e9;       % [Pa]
E_2024 = 73e9;       % [Pa]
E_6061 = 69e9;       % [Pa]

Re_7075 = 503e6;     % [Pa] yield
Re_2024 = 324e6;     % [Pa] yield
Re_6061 = 275e6;     % [Pa] yield

SF = 1.5;            % safety factor

%   bottom skin, bottom cap, web, top cap, top skin
E = [E_7075, E_2024, E_6061, E_2024, E_7075];  % [Pa]

sigma_allow = [Re_7075, Re_2024, Re_6061, Re_2024, Re_7075] / SF;  % [Pa]

%% ------------------------------------------------------------------------
%  INTERNAL FORCES IN CRITICAL SECTION
% -------------------------------------------------------------------------
N =  0;     % [N] axial (positive tension)
M = 219.1779307e3;       % [Nm] bending moment

%% ------------------------------------------------------------------------
%  BASIC SECTION PROPERTIES (AREAS & "B" STIFFNESS FOR N)
% -------------------------------------------------------------------------
B_i = A .* E;          % [N]   (Bi = Ai*Ei)
B   = sum(B_i);        % [N]

% total height and centroid positions of each layer (from bottom)
H      = sum(t);
y_cent = zeros(1,nLayers);

y_running = 0;
for i = 1:nLayers
    y_cent(i) = y_running + t(i)/2;
    y_running = y_running + t(i);
end

%% ------------------------------------------------------------------------
%  NEUTRAL AXIS POSITION (bending) – weighted by Bi
% -------------------------------------------------------------------------
y_NA = sum(B_i .* y_cent) / B   % [m]

term1 = sum(B_i .* t);
term2 = sum(B_i(2:end) .* cumsum(t(1:end-1)));
y_NA_check = (term1 + 2*term2) / (2*sum(B_i))

% distance of each layer centroid to NA (positive upwards)
y_iN = y_cent - y_NA;            % [m]

%% ------------------------------------------------------------------------
%  SECOND MOMENT OF AREA & BENDING STIFFNESS D
% -------------------------------------------------------------------------
I_cent = (b .* t.^3) / 12;       % I of each layer about its own centroid
I_i    = I_cent + A .* (y_iN.^2);% parallel axis to NA

D_i = E .* I_i;                  % [Nm^2]
D   = sum(D_i);                  % [Nm^2]  (overall bending stiffness)

%% ------------------------------------------------------------------------
%  NORMAL STRESSES  σ_i = E_i*( N/B + M*y_i / D )
% -------------------------------------------------------------------------
epsilon_N = N / B;               % common axial strain
sigma_i   = E .* ( epsilon_N + (M .* y_iN) / D );  % [Pa] at layer centroids

% extreme fibres (bottom & top of whole section)
y_top = H - y_NA;
y_bot = -y_NA;

sigma_bot = E(1)   * (epsilon_N + M*y_bot/D);   % bottom skin material
sigma_top = E(end) * (epsilon_N + M*y_top/D);   % top skin material

%% ------------------------------------------------------------------------
%  UTILISATION (NORMAL STRESSES ONLY)
% -------------------------------------------------------------------------
util_sigma_layers = abs(sigma_i)   ./ sigma_allow;
util_sigma_bot    = abs(sigma_bot) / sigma_allow(1);
util_sigma_top    = abs(sigma_top) / sigma_allow(end);

%% ------------------------------------------------------------------------
%  PRINT RESULTS
% -------------------------------------------------------------------------
fprintf('--- MULTILAYER SECTION STRENGTH CHECK (NORMAL STRESS) ---\n');
fprintf('Neutral axis from bottom y_NA = %.4f m\n', y_NA);
fprintf('Axial stiffness B        = %.3e N\n', B);
fprintf('Bending stiffness D      = %.3e Nm^2\n\n', D);

for i = 1:nLayers
    fprintf(['Layer %d: A = %.3e m^2, y_cent = %.4f m, ' ...
             'sigma = %+8.2f MPa, util = %.2f\n'], ...
        i, A(i), y_cent(i), sigma_i(i)/1e6, util_sigma_layers(i));
end

fprintf('\nExtreme fibres:\n');
fprintf('Bottom sigma = %+8.2f MPa, utilisation = %.2f\n', ...
    sigma_bot/1e6, util_sigma_bot);
fprintf('Top    sigma = %+8.2f MPa, utilisation = %.2f\n', ...
    sigma_top/1e6, util_sigma_top);

%% -----------------------------------------------------------
%  NORMAL STRESS DISTRIBUTION PLOT σ(y)
% ------------------------------------------------------------
% layer interfaces from bottom: 0, t1, t1+t2, ...
y_iface = [0, cumsum(t)];

% global axial strain (already epsilon_N)
nPts = 300;                        % resolution of the graph
y_vec = linspace(0, H, nPts);      % from bottom (0) to top (H)
sigma_y = zeros(size(y_vec));      % [Pa]

for k = 1:nPts
    y = y_vec(k);

    % find which layer this y belongs to
    idx = find(y >= y_iface(1:end-1) & y <= y_iface(2:end), 1, 'first');

    % distance from NA
    y_rel = y - y_NA;

    % strain and stress in that layer
    epsilon_y = epsilon_N + (M * y_rel) / D;
    sigma_y(k) = E(idx) * epsilon_y;
end

% convert to MPa and flip sign (for aircraft convention: tension at bottom)
sigma_MPa = -sigma_y / 1e6;

figure; hold on; grid on;
plot(sigma_MPa, y_vec, 'LineWidth', 1.8);

xline(0, '--k', 'LineWidth', 1);   % σ = 0
yline(y_NA, ':k', 'LineWidth', 1); % NA position

xlabel('\sigma [MPa]');
ylabel('y [m] (from bottom of spar)');
title('Normal stress distribution \sigma(y) in spar cross-section');

set(gca, 'YDir', 'normal');  % bottom at y=0, top at y=H

