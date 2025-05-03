clc;
clear;
close;

%% informazioni sulla trave
L = 10;        % Lunghezza della trave (m)
E = 210e9;     % Modulo di Young (Pa)
I = 8.33e-6;   % Momento d'inerzia della sezione (m^4)
rho = 7800;    % Densità del materiale (kg/m^3)
A = 0.01;      % Area della sezione trasversale (m^2)

N=14; %% N elementi, N+1 nodi

L_el=L/N; 

% Matrice di rigidezza per un elemento 1D (flessione), 2 GDL per nodo, per ogni elemento 4 GDL
k_el = (E * I / L_el^3) * [12, 6*L_el, -12, 6*L_el; 
                      6*L_el, 4*L_el^2, -6*L_el, 2*L_el^2;
                      -12, -6*L_el, 12, -6*L_el;
                      6*L_el, 2*L_el^2, -6*L_el, 4*L_el^2];

% Matrice di massa per un elemento 1D (flessione)
m_el = (rho * A * L_el / 420) * [156, 22*L_el, 54, -13*L_el; 
                            22*L_el, 4*L_el^2, 13*L_el, -3*L_el^2;
                            54, 13*L_el, 156, -22*L_el;
                            -13*L_el, -3*L_el^2, -22*L_el, 4*L_el^2];

% Assemblaggio delle matrici di rigidezza e massa globali
K = zeros(2*(N+1), 2*(N+1));  % Matrice di rigidezza globale, 2(N+1) GDL
M = zeros(2*(N+1), 2*(N+1));  % Matrice di massa globale

for p = 1:N
    nodes = [2*p-1, 2*p, 2*p+1, 2*p+2]; 
    K(nodes, nodes) = K(nodes, nodes) + k_el;
    M(nodes, nodes) = M(nodes, nodes) + m_el;
end

% Condizioni al contorno: la trave è appoggiata su entrambi i lati (vincolo traslazione)
% La trave è appoggiata sui nodi 1 e 2n+1 (posizioni in x)

% Condizioni al contorno (nodi fissi)
fixed_nodes = [1, 2*N+1];   % Nodo 1 e nodo 2 sono vincolati alla traslazione
free_nodes = setdiff(1:2*(N+1), fixed_nodes);  % Nodi liberi

% Ridurre il sistema per i nodi liberi
K = K(free_nodes, free_nodes);    % Matrice di rigidezza per i nodi liberi
M = M(free_nodes, free_nodes);    % Matrice di massa per i nodi liberi

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%        PROPORTIONAL STRUCTURAL DAMPING               %%%%%%%%%%%%%
beta=0.05;
D_1 = beta*K;
[PHI_PROP_DUMP, LAMBDA_PROP_DUMP] = eig(K+1j*D_1, M);
[LAMBDA_PROP_DUMP,ind_2]=sort(diag(LAMBDA_PROP_DUMP));
LAMBDA_PROP_DUMP=diag(LAMBDA_PROP_DUMP);
PHI_PROP_DUMP=PHI_PROP_DUMP(:,ind_2);

for i = 1:size(PHI_PROP_DUMP, 2)  % Itera su ogni autovettore
    v = PHI_PROP_DUMP(:, i);  % Estrai l'autovettore i-esimo
    norm_factor = sqrt(v' * M * v);  % Calcola il fattore di normalizzazione
    PHI_PROP_DUMP(:, i) = v / norm_factor;  % Normalizza l'autovettore
end

%%
%%%%%%%%%%%%              FRF                                    %%%%%%%%%%
L = 50000;
omega_max = sqrt(real(LAMBDA_PROP_DUMP(10,10)));
% Range di frequenze per il calcolo della FRF
omega = linspace(0, 1.3*omega_max, L);

% Numero di modi che considero nella sommatoria 
n = N/2; 

% Risposta nel nodo j
j = 3;

% Forza applicata nel nodo k
k = 3; 

%FRF
FRF_PROP_DUMP = zeros(1,L);
for r=1:n
  FRF_PROP_DUMP = FRF_PROP_DUMP+...
      ((PHI_PROP_DUMP(j,r)*PHI_PROP_DUMP(k,r))./...
        (LAMBDA_PROP_DUMP(r,r)-(omega.^2)));
end


%% grafica FRF


figure;
% 
% semilogy(omega, abs(FRF_PROP_DUMP), 'r', 'LineWidth', 1.2);
% xlabel('\omega [rad/s]');
% ylabel('|\alpha_{j,k}| [m/N]');
% title('FRF in scala log-log');
% grid on;


% Plot della FRF  
plot(omega, abs(FRF_PROP_DUMP), 'r--', 'LineWidth', 1.2);
hold on;

% Aggiunta di titoli e etichette
title('Recettanza \alpha_{3,3} per trave');
xlabel('\omega [rad/s]');
ylabel('|\alpha_{3,3}| [m/N]');
% Aggiunta griglia e leggenda
grid on;


omega_nat_vec = sqrt(diag(real(LAMBDA_PROP_DUMP(1:3,1:3))))
A1_3_3= PHI_PROP_DUMP(3,1)*PHI_PROP_DUMP(3,1)
A2_3_3= PHI_PROP_DUMP(3,2)*PHI_PROP_DUMP(3,2)
A3_3_3= PHI_PROP_DUMP(3,3)*PHI_PROP_DUMP(3,3)

%%
% FRF completa già calcolata: omega, FRF_PROP_DUMP

% 1. Trova il massimo (prima risonanza per semplicità)
% [~, peak_idx] = max(abs(FRF_PROP_DUMP));
% w_res = omega(peak_idx);  % Frequenza di risonanza stimata

[pks, locs] = findpeaks_basic(abs(FRF_PROP_DUMP));

% Seleziona il secondo picco risonanza
peak_idx = locs(3);
w_res = omega(peak_idx);  % Frequenza di risonanza stimata


% 2. Seleziona zona intorno alla risonanza
window_pct = 0.1; % 1% attorno al picco
intorno = find(omega > (1-window_pct)*w_res & omega < (1+window_pct)*w_res);

omega_fit = omega(intorno);
FRF_fit = FRF_PROP_DUMP(intorno);
FRF_fit = FRF_fit(:);
% 3. Circle fit nel piano complesso
x = real(FRF_fit);
y = imag(FRF_fit);

A = [-2*x, -2*y, ones(length(x),1)];
b = -(x.^2 + y.^2);
params = A\b;

xc = params(1); yc = params(2);
R = sqrt(xc^2 + yc^2 - params(3));
D = 2 * R;

% === Interpolazione spline complessa della FRF ===
FRF_spline = @(w) interp1(omega_fit, FRF_fit, w, 'spline');

% === Metodo 1: Max parte immaginaria ===
Im_FRF_spline = @(w) imag(FRF_spline(w));
omega_n_1 = fminbnd(@(w) Im_FRF_spline(w), omega_fit(1), omega_fit(end));
z1 = FRF_spline(omega_n_1);

% === Metodo 2: Punto con Re = 0 ===
Re_FRF_spline = @(w) real(FRF_spline(w));
omega_n_2 = fzero(Re_FRF_spline, w_res);  % Guess vicino al picco
z2 = FRF_spline(omega_n_2);

% === Metodo 3: Distanza massima dall’origine (modulo massimo) ===
Mod_FRF_spline = @(w) abs(FRF_spline(w));
omega_n_3 = fminbnd(@(w) -Mod_FRF_spline(w), omega_fit(1), omega_fit(end));
z3 = FRF_spline(omega_n_3);

% === Plot della Circle Fit con i 3 punti stimati ===
theta = linspace(0, 2*pi, 300);
circle = xc + 1i*yc + R * exp(1i*theta);

figure;
plot(real(FRF_fit), imag(FRF_fit), 'bo', 'DisplayName', 'FRF zona risonanza'); hold on;
plot(real(circle), imag(circle), 'k-', 'LineWidth', 1.5, 'DisplayName', 'Circle fit');
plot(real(z1), imag(z1), 'r+', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', '\omega_1 max Im');
plot(real(z2), imag(z2), 'b+', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', '\omega_1 Re = 0');
plot(real(z3), imag(z3), 'm+', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', '\omega_1 max |FRF|');
plot(xc, yc, 'g*', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Centro cerchio');

axis equal;
xlabel('Re\{\alpha\}');
ylabel('Im\{\alpha\}');
title('Stima della frequenza naturale \omega_1 con 3 metodi (circle fit)');
legend show;
grid on;

% === Stampa dei risultati in formato LaTeX ===
fprintf('\n=== Frequenze naturali stimate sulla circle fit ===\n');
fprintf('Metodo 1 (max Im):       ω₃ = %.4f rad/s\n', omega_n_1);
fprintf('Metodo 2 (Re = 0):       ω₃  = %.4f rad/s\n', omega_n_2);
fprintf('Metodo 3 (max |FRF|):    ω₃  = %.4f rad/s\n', omega_n_3);

% Shift della FRF rispetto al centro del cerchio
C = xc + 1i*yc;
FRF_shifted = FRF_fit - C;

% Angoli rispetto al centro (unwrap per continuità)
theta_FRF = unwrap(angle(FRF_shifted));

% Interpolazione per theta = 0 (destra del cerchio)
omega_a = interp1(theta_FRF, omega_fit, 0, 'linear');

% Interpolazione per theta = pi (sinistra del cerchio)
omega_b = interp1(theta_FRF, omega_fit, -pi, 'linear');

% Output
fprintf('\nomega_a (θ = 0):    %.4f rad/s  [destra del cerchio]\n', omega_a);
fprintf('omega_b (θ = π):    %.4f rad/s  [sinistra del cerchio]\n', omega_b);


% === Calcolo dello smorzamento ===
zeta = (omega_a - omega_b) / omega_n_1;  % Rapporto delle frequenze

% === Visualizzazione dei punti omega_a e omega_b sulla circonferenza ===
z_a = FRF_spline(omega_a); % Punto omega_a
z_b = FRF_spline(omega_b); % Punto omega_b

% Plot aggiornata con i punti omega_a e omega_b
figure;
plot(real(FRF_fit), imag(FRF_fit), 'bo', 'DisplayName', 'FRF zona risonanza'); hold on;
plot(real(circle), imag(circle), 'k-', 'LineWidth', 1.5, 'DisplayName', 'Circle fit');
plot(real(z1), imag(z1), 'r+', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', '\omega_n max Im');
plot(real(z2), imag(z2), 'b+', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', '\omega_n Re = 0');
plot(real(z3), imag(z3), 'm+', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', '\omega_n max |FRF|');
plot(real(z_a), imag(z_a), 'g*', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', '\omega_a');
plot(real(z_b), imag(z_b), 'g*', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', '\omega_b');
plot(xc, yc, 'g+', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Centro cerchio');

axis equal;
xlabel('Re\{\alpha\}');
ylabel('Im\{\alpha\}');
title('Stima dello smorzamento modale ');
legend show;
grid on;


fprintf('Smorzamento stimato:    η₃ = %.4f\n', zeta);



%% Calcolo costante modale

A33_1 = 2*R*(omega_n_1^2)*zeta;
fprintf('\nDiametro:   D = %f\n', D);
fprintf('Costante modale stimata:    ₃A₃,₃ = %f\n', A33_1);




%%%
z_max = z3;

%% Plot del Circle Fit con diametro principale

theta = linspace(0, 2*pi, 200);
circle = xc + 1i*yc + R * exp(1i*theta);

% Punto diametralmente opposto a z_max
z_opposto = 2*(xc + 1i*yc) - z_max;

figure;
hold on; grid on; axis equal;

% Plot FRF nel piano complesso (zona di risonanza)
plot(real(FRF_fit), imag(FRF_fit), 'bo', 'DisplayName', 'FRF zona risonanza');

% Circle fit
plot(real(circle), imag(circle), 'k-', 'LineWidth', 1.5, 'DisplayName', 'Circle fit');

% Centro del cerchio
plot(xc, yc, 'g+', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Centro cerchio');

% Punto di risonanza (massimo |FRF|)
plot(real(z3), imag(z3), 'm+', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', '\omega_n max |FRF|');

% Punto max Im
plot(real(z1), imag(z1), 'r+', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', '\omega_n max Im');

% Punto Re = 0
plot(real(z2), imag(z2), 'b+', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', '\omega_n Re = 0');

% Punti omega_a e omega_b
plot(real(z_a), imag(z_a), 'g*', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', '\omega_a (\theta = 0)');
plot(real(z_b), imag(z_b), 'c*', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', '\omega_b (\theta = \pi)');

% Diametro principale
z_opposto = 2*(xc + 1i*yc) - z3;
plot([real(z3), real(z_opposto)], [imag(z3), imag(z_opposto)], ...
     'm--', 'LineWidth', 2, 'DisplayName', 'Diametro principale');

xlabel('Re\{\alpha\}');
ylabel('Im\{\alpha\}');
title('Circle Fit FRF: stimatori \omega_3, \eta_3, {}_{3}A_{3,3}');
legend('Location', 'bestoutside');



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


