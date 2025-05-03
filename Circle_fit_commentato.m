clc;
clear;
close all;

%% ========================================================================
%                    DATI DELLA TRAVE E MESH
%==========================================================================
L = 10;          % Lunghezza trave [m]
E = 210e9;       % Modulo di Young [Pa]
I = 8.33e-6;     % Momento d'inerzia [m^4]
rho = 7800;      % Densità [kg/m^3]
A = 0.01;        % Area sezione trasversale [m^2]

N = 14;          % Numero di elementi (N+1 nodi)
L_el = L / N;    % Lunghezza di ciascun elemento

%% ========================================================================
%                   MATRICI ELEMENTARI (rigidezza e massa)
%==========================================================================
k_el = (E * I / L_el^3) * [12, 6*L_el, -12, 6*L_el; 
                           6*L_el, 4*L_el^2, -6*L_el, 2*L_el^2;
                          -12, -6*L_el, 12, -6*L_el;
                           6*L_el, 2*L_el^2, -6*L_el, 4*L_el^2];

m_el = (rho * A * L_el / 420) * [156, 22*L_el, 54, -13*L_el; 
                                 22*L_el, 4*L_el^2, 13*L_el, -3*L_el^2;
                                 54, 13*L_el, 156, -22*L_el;
                                -13*L_el, -3*L_el^2, -22*L_el, 4*L_el^2];

%% ========================================================================
%                   ASSEMBLAGGIO GLOBALI
%==========================================================================
K = zeros(2*(N+1));
M = zeros(2*(N+1));

for p = 1:N
    nodes = [2*p-1, 2*p, 2*p+1, 2*p+2]; 
    K(nodes, nodes) = K(nodes, nodes) + k_el;
    M(nodes, nodes) = M(nodes, nodes) + m_el;
end

%% ========================================================================
%                   VINCOLI (condizioni al contorno)
%==========================================================================
fixed_nodes = [1, 2*N+1];  % Estremità vincolate
free_nodes = setdiff(1:2*(N+1), fixed_nodes);

K = K(free_nodes, free_nodes);
M = M(free_nodes, free_nodes);

%% ========================================================================
%                   SMORZAMENTO STRUTTURALE PROPORZIONALE
%==========================================================================
beta = 0.05;
D_1 = beta * K;

[PHI_PROP_DUMP, LAMBDA_PROP_DUMP] = eig(K + 1j*D_1, M);
[LAMBDA_PROP_DUMP, ind] = sort(diag(LAMBDA_PROP_DUMP));
LAMBDA_PROP_DUMP = diag(LAMBDA_PROP_DUMP);
PHI_PROP_DUMP = PHI_PROP_DUMP(:, ind);

% Normalizzazione modale (massa unitaria)
for i = 1:size(PHI_PROP_DUMP, 2)
    v = PHI_PROP_DUMP(:, i);
    PHI_PROP_DUMP(:, i) = v / sqrt(v' * M * v);
end

%% ========================================================================
%                   FRF (Funzione di Risposta in Frequenza)
%==========================================================================
L = 50000;
omega_max = sqrt(real(LAMBDA_PROP_DUMP(10,10)));
omega = linspace(0, 1.3*omega_max, L);

n = N / 2;
j = 3;  % Nodo di risposta
k = 3;  % Nodo di applicazione forza

FRF_PROP_DUMP = zeros(1, L);
for r = 1:n
    FRF_PROP_DUMP = FRF_PROP_DUMP + ...
        (PHI_PROP_DUMP(j, r) * PHI_PROP_DUMP(k, r)) ./ ...
        (LAMBDA_PROP_DUMP(r, r) - omega.^2);
end

%% ========================================================================
%                   PLOT DELLA FRF
%==========================================================================
figure;
plot(omega, abs(FRF_PROP_DUMP), 'r--', 'LineWidth', 1.2);
title('Recettanza \alpha_{3,3} per trave');
xlabel('\omega [rad/s]');
ylabel('|\alpha_{3,3}| [m/N]');
grid on;

omega_nat_vec = sqrt(diag(real(LAMBDA_PROP_DUMP(1:3, 1:3))));
A1_3_3 = PHI_PROP_DUMP(3,1)^2;
A2_3_3 = PHI_PROP_DUMP(3,2)^2;
A3_3_3 = PHI_PROP_DUMP(3,3)^2;

%% ========================================================================
%                   STIMA FREQUENZA NATURALE (CIRCLE FIT)
%==========================================================================
[pks, locs] = findpeaks_basic(abs(FRF_PROP_DUMP));
peak_idx = locs(3);  % Terza risonanza
w_res = omega(peak_idx);

% Intervallo per interpolazione
window_pct = 0.1;
intorno = find(omega > (1-window_pct)*w_res & omega < (1+window_pct)*w_res);
omega_fit = omega(intorno);
FRF_fit = FRF_PROP_DUMP(intorno);

% Circle fit nel piano complesso
x = real(FRF_fit);
x = x(:);
y = imag(FRF_fit);
y = y(:);
A = [-2*x, -2*y, ones(length(x),1)];
b = -(x.^2 + y.^2);
params = A\b;

xc = params(1); yc = params(2);
R = sqrt(xc^2 + yc^2 - params(3));
D = 2 * R;
C = xc + 1i*yc;

% Interpolazione spline complessa
FRF_spline = @(w) interp1(omega_fit, FRF_fit, w, 'spline');

omega_n_1 = fminbnd(@(w) imag(FRF_spline(w)), omega_fit(1), omega_fit(end));
omega_n_2 = fzero(@(w) real(FRF_spline(w)), w_res);
omega_n_3 = fminbnd(@(w) -abs(FRF_spline(w)), omega_fit(1), omega_fit(end));

z1 = FRF_spline(omega_n_1);
z2 = FRF_spline(omega_n_2);
z3 = FRF_spline(omega_n_3);

% Stampa risultati frequenze naturali
fprintf('\n=== Frequenze naturali stimate sulla circle fit ===\n');
fprintf('Metodo 1 (max Im):       ω₃ = %.4f rad/s\n', omega_n_1);
fprintf('Metodo 2 (Re = 0):       ω₃ = %.4f rad/s\n', omega_n_2);
fprintf('Metodo 3 (max |FRF|):    ω₃ = %.4f rad/s\n', omega_n_3);

%% ========================================================================
%                   STIMA SMORZAMENTO MODALE
%==========================================================================
FRF_shifted = FRF_fit - C;
theta_FRF = unwrap(angle(FRF_shifted));

omega_a = interp1(theta_FRF, omega_fit, 0, 'linear');
omega_b = interp1(theta_FRF, omega_fit, -pi, 'linear');

zeta = (omega_a - omega_b) / omega_n_1;

fprintf('\nomega_a (θ = 0):    %.4f rad/s\n', omega_a);
fprintf('omega_b (θ = π):    %.4f rad/s\n', omega_b);
fprintf('Smorzamento stimato:    η₃ = %.4f\n', zeta);

% Costante modale stimata
A33_1 = 2 * R * (omega_n_1^2) * zeta;
fprintf('\nDiametro:   D = %f\n', D);
fprintf('Costante modale stimata:    ₃A₃,₃ = %f\n', A33_1);

%% ========================================================================
%                   PLOT FINALE: CIRCLE FIT COMPLETO
%==========================================================================
z_a = FRF_spline(omega_a);
z_b = FRF_spline(omega_b);
z_max = z3;
z_opposto = 2*C - z_max;

theta = linspace(0, 2*pi, 200);
circle = xc + 1i*yc + R * exp(1i*theta);

figure;
hold on; grid on; axis equal;

plot(real(FRF_fit), imag(FRF_fit), 'bo', 'DisplayName', 'FRF zona risonanza');
plot(real(circle), imag(circle), 'k-', 'LineWidth', 1.5, 'DisplayName', 'Circle fit');
plot(xc, yc, 'g+', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Centro cerchio');

plot(real(z3), imag(z3), 'm+', 'LineWidth', 2, 'DisplayName', '\omega_n max |FRF|');
plot(real(z1), imag(z1), 'r+', 'LineWidth', 2, 'DisplayName', '\omega_n max Im');
plot(real(z2), imag(z2), 'b+', 'LineWidth', 2, 'DisplayName', '\omega_n Re = 0');
plot(real(z_a), imag(z_a), 'g*', 'LineWidth', 2, 'DisplayName', '\omega_a (\theta = 0)');
plot(real(z_b), imag(z_b), 'c*', 'LineWidth', 2, 'DisplayName', '\omega_b (\theta = \pi)');
plot([real(z3), real(z_opposto)], [imag(z3), imag(z_opposto)], ...
     'm--', 'LineWidth', 2, 'DisplayName', 'Diametro principale');

xlabel('Re\{\alpha\}');
ylabel('Im\{\alpha\}');
title('Circle Fit FRF: stimatori \omega_3, \eta_3, {}_{3}A_{3,3}');
legend('Location', 'bestoutside');

%% ========================================================================
%                   FUNZIONE UTILE: FINDPEAKS BASIC
%==========================================================================
function [pks, locs] = findpeaks_basic(y)
    pks = [];
    locs = [];
    for i = 2:length(y)-1
        if y(i) > y(i-1) && y(i) > y(i+1)
            pks(end+1) = y(i);
            locs(end+1) = i;
        end
    end
end
