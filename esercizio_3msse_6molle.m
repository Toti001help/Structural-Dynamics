clc;
clear;
close all;
N = 3; % GDL
m1=1; % Kg
m2=0.95; % Kg
m3=1.05; % Kg

k=1.0e3; % N/m

M=[m1, 0, 0;      % Mass Matrix
    0, m2, 0;
    0, 0, m3];

K=[k+k+k, -k, -k; % Stiffness Matrix
    -k, k+k+k, -k;
    -k, -k, k+k+k];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%     UNDAMPED     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[PHI_UN, LAMBDA_UN] = eig(K, M);
[LAMBDA_UN,ind_1]=sort(diag(LAMBDA_UN));
LAMBDA_UN=diag(LAMBDA_UN);
PHI_UN=PHI_UN(:,ind_1);

for i = 1:size(PHI_UN, 2)  % Itera su ogni autovettore
    v = PHI_UN(:, i);  % Estrai l'autovettore i-esimo
    norm_factor = sqrt(v' * M * v);  % Calcola il fattore di normalizzazione
    PHI_UN(:, i) = v / norm_factor;  % Normalizza l'autovettore
end

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%        NON PROPORTIONAL STRUCTURAL DAMPING               %%%%%%%%%
d=0.3*k;
D_2 = zeros(3,3);
D_2(1,1) = d;
[PHI_NOPROP, LAMBDA_NO_PROP] = eig(K+(1j*D_2), M);
[LAMBDA_NO_PROP,ind_3]=sort(diag(LAMBDA_NO_PROP));
LAMBDA_NO_PROP=diag(LAMBDA_NO_PROP);
PHI_NOPROP=PHI_NOPROP(:,ind_3);

for i = 1:size(PHI_NOPROP, 2)  % Itera su ogni autovettore complesso
    v = PHI_NOPROP(:, i);  % Estrai l'autovettore i-esimo
    norm_factor = sqrt(conj(v)' * M * v);  % Normalizzazione Hermitiana
    PHI_NOPROP(:, i) = v / norm_factor;  % Normalizza l'autovettore
end

%%                                                                       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%              FRF                                    %%%%%%%%%%
L = 1000;
omega_max = sqrt(LAMBDA_UN(3,3));
% Range di frequenze per il calcolo della FRF
omega = linspace(0, 1.3*omega_max, L);

% Risposta nel nodo j
j = 1;

% Forza applicata nel nodo i
k = 1; 

% Calcolo della FRF per i tre modelli: senza smorzamento, proporzionale, 
% non proporzionale
%inizializzo risposte 
FRF_UN = zeros(1,L);
FRF_PROP_DUMP = zeros(1,L);
FRF_NO_PROP = zeros(1,L);

FRF_UN_rip = zeros(1,L);
FRF_PROP_DUMP_rip = zeros(1,L);
FRF_NO_PROP_rip = zeros(1,L);

% Caso senza smorzamento
%D_UN = zeros(3, 3);  % Matrice di smorzamento nulla

% for i = 1:L
%     %% Caso senza smorzamento
%     H1 = inv(K+1j*D_UN-((omega(i))^2)*M);
%     FRF_UN(i) = H1(j,k);
% 
% 
%     %% Caso smorzamento Proporzionale
%     H2 = inv(K+1j*D_1-((omega(i))^2)*M);
%     FRF_PROP_DUMP(i) = H2(j,k);
% 
% 
%     %% Caso smorzamento non proporzionale
%     H3 = inv(K+1j*D_2-((omega(i))^2)*M);
%     FRF_NO_PROP(i) = H3(j,k);
% end

for r=1:N
%% Caso senza smorzamento
  FRF_UN_rip = FRF_UN_rip+((PHI_UN(j,r)*PHI_UN(k,r))./...
        (LAMBDA_UN(r,r)-(omega.^2)));
%% Caso smorzamento Proporzionale
  FRF_PROP_DUMP_rip = FRF_PROP_DUMP_rip+...
      ((PHI_PROP_DUMP(j,r)*PHI_PROP_DUMP(k,r))./...
        (LAMBDA_PROP_DUMP(r,r)-(omega.^2)));
%% Caso smorzamento non proporzionale
  FRF_NO_PROP_rip = FRF_NO_PROP_rip+((PHI_NOPROP(j,r)*PHI_NOPROP(k,r))./...
        (LAMBDA_NO_PROP(r,r)-(omega.^2)));
end





% 
% L = 10000;
% omega_max = sqrt(LAMBDA_UN(3,3));
% % Range di frequenze per il calcolo della FRF
% omega2 = linspace(0, 1.3*omega_max, L);
% 
% % Risposta nel nodo j
% j = 1;
% 
% % Forza applicata nel nodo i
% k = 1; 
% 
% % Calcolo della FRF per i tre modelli: senza smorzamento, proporzionale, 
% % non proporzionale
% %inizializzo risposte 
% FRF_UN2 = zeros(1,L);
% FRF_PROP_DUMP2 = zeros(1,L);
% FRF_NO_PROP2 = zeros(1,L);
% 
% FRF_UN_rip2 = zeros(1,L);
% FRF_PROP_DUMP_rip2 = zeros(1,L);
% FRF_NO_PROP_rip2 = zeros(1,L);

% Caso senza smorzamento
%D_UN = zeros(3, 3);  % Matrice di smorzamento nulla


% for r=1:N
% %% Caso senza smorzamento
%   FRF_UN_rip2 = FRF_UN_rip2+((PHI_UN(j,r)*PHI_UN(k,r))./...
%         (LAMBDA_UN(r,r)-(omega2.^2)));
% %% Caso smorzamento Proporzionale
%   FRF_PROP_DUMP_rip2 = FRF_PROP_DUMP_rip2+...
%       ((PHI_PROP_DUMP(j,r)*PHI_PROP_DUMP(k,r))./...
%         (LAMBDA_PROP_DUMP(r,r)-(omega2.^2)));
% %% Caso smorzamento non proporzionale
%   FRF_NO_PROP_rip2 = FRF_NO_PROP_rip2+((PHI_NOPROP(j,r)*PHI_NOPROP(k,r))./...
%         (LAMBDA_NO_PROP(r,r)-(omega2.^2)));
% end





% 
% L = 100000;
% omega_max = sqrt(LAMBDA_UN(3,3));
% % Range di frequenze per il calcolo della FRF
% omega3 = linspace(0, 1.3*omega_max, L);
% 
% % Risposta nel nodo j
% j = 1;
% 
% % Forza applicata nel nodo i
% k = 1; 
% 
% % Calcolo della FRF per i tre modelli: senza smorzamento, proporzionale, 
% % non proporzionale
% %inizializzo risposte 
% FRF_UN = zeros(1,L);
% FRF_PROP_DUMP = zeros(1,L);
% FRF_NO_PROP = zeros(1,L);
% 
% FRF_UN_rip3 = zeros(1,L);
% FRF_PROP_DUMP_rip3 = zeros(1,L);
% FRF_NO_PROP_rip3 = zeros(1,L);
% 
% for r=1:N
% %% Caso senza smorzamento
%   FRF_UN_rip3 = FRF_UN_rip3+((PHI_UN(j,r)*PHI_UN(k,r))./...
%         (LAMBDA_UN(r,r)-(omega3.^2)));
% %% Caso smorzamento Proporzionale
%   FRF_PROP_DUMP_rip3 = FRF_PROP_DUMP_rip3+...
%       ((PHI_PROP_DUMP(j,r)*PHI_PROP_DUMP(k,r))./...
%         (LAMBDA_PROP_DUMP(r,r)-(omega3.^2)));
% %% Caso smorzamento non proporzionale
%   FRF_NO_PROP_rip3 = FRF_NO_PROP_rip3+((PHI_NOPROP(j,r)*PHI_NOPROP(k,r))./...
%         (LAMBDA_NO_PROP(r,r)-(omega3.^2)));
% end
% 
% 
% 
% 






 % figure;
 % plot(omega3, abs(FRF_UN_rip3), 'b-', 'LineWidth', 2);
 % hold on
 % plot(omega2, abs(FRF_UN_rip2), 'r-', 'LineWidth', 2);
 % 
 % %plot(omega, abs(FRF_UN_rip2), 'g--', 'LineWidth', 2);
 % 
 % title('FRF senza Smorzamento, diverso passo in frequenza');
 % xlabel('\omega [rad/s]');
 % ylabel('\alpha{1,1}');
 % grid on;
 % 
 % legend( 'N_{\omega}=10000',...
 %     'N_{\omega}=100000', ...
 %    'Location', 'Best');
 % xlim([min(omega), max(omega)]);
 % ylim( [0, 3*max(FRF_UN_rip2)])


figure;

% Plot della FRF senza smorzamento (linea rossa tratteggiata sottile)
plot(omega, abs(FRF_UN_rip), 'r--', 'LineWidth', 1.2);
hold on;

% Plot della FRF con smorzamento proporzionale (linea verde continua)
plot(omega, abs(FRF_PROP_DUMP_rip), 'g-', 'LineWidth', 2);

% Plot della FRF con smorzamento non proporzionale (linea blu continua)
plot(omega, abs(FRF_NO_PROP_rip), 'b-', 'LineWidth', 2);

% Aggiunta di titoli e etichette
title('Funzione di cross recettanza \alpha_{1,1}');
xlabel('\omega [rad/s]');
ylabel('|\alpha_{1,1}| [m/N]');

% Aggiunta griglia e leggenda
grid on;
legend('Senza Smorzamento', 'Smorzamento Proporzionale', 'Smorzamento Non Proporzionale', ...
    'Location', 'Best');

% Miglioramento grafico
xlim([min(omega), max(omega)]);
ylim( [0, 1.3*max(max(abs(FRF_PROP_DUMP_rip)),max(abs(FRF_NO_PROP_rip)))])













% FRF completa già calcolata: omega, FRF_PROP_DUMP

% 1. Trova il massimo (prima risonanza per semplicità)
% [~, peak_idx] = max(abs(FRF_NO_PROP_rip));
% w_res = omega(peak_idx);  % Frequenza di risonanza stimata

[pks, locs] = findpeaks_basic(abs(FRF_PROP_DUMP_rip));

% Seleziona il secondo picco risonanza
peak_idx = locs(2);
w_res = omega(peak_idx);  % Frequenza di risonanza stimata


% 2. Seleziona zona intorno alla risonanza
window_pct = 0.1; % 1% attorno al picco
intorno = find(omega > (1-window_pct)*w_res & omega < (1+window_pct)*w_res);

omega_fit = omega(intorno);
FRF_fit = FRF_PROP_DUMP_rip(intorno);
FRF_fit = FRF_fit(:);
% 3. Circle fit nel piano complesso
x = real(FRF_fit);
y = imag(FRF_fit);

A = [-2*x, -2*y, ones(length(x),1)];
b = -(x.^2 + y.^2);
params = A\b;

xc = params(1); yc = params(2);
R = sqrt(xc^2 + yc^2 - params(3));

% 4. Calcolo smorzamento ζ
% Punto massimo (vicino al bordo del cerchio)
z_max = FRF_PROP_DUMP_rip(peak_idx);
dist_to_center = abs(z_max - (xc + 1i*yc));
zeta = dist_to_center / R;

% 5. Ampiezza modale ≈ raggio
A_modale = R;

% 6. Output dei risultati
fprintf('\nFrequenza di risonanza stimata: %.2f rad/s\n', w_res);
fprintf('Smorzamento stimato (zeta): %.4f\n', zeta);
fprintf('Ampiezza modale (raggio del cerchio): %.3e [unità FRF]\n', A_modale);

% 7. Grafico circle fit
theta = linspace(0, 2*pi, 200);
circle = xc + 1i*yc + R * exp(1i*theta);

figure;
plot(real(FRF_fit), imag(FRF_fit), 'bo'); hold on;
plot(real(circle), imag(circle), 'r-', 'LineWidth', 1.5);
plot(real(z_max), imag(z_max), 'kx', 'LineWidth', 2, 'MarkerSize', 10);
plot(xc, yc, 'g+', 'MarkerSize', 10, 'LineWidth', 2);
axis equal;
xlabel('Re'); ylabel('Im');
title('Circle Fit + Stima ζ, ω_n, A');
legend('FRF (zona risonanza)', 'Cerchio', 'Picco FRF', 'Centro cerchio');
grid on;



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





% figure;
% plot(omega, abs(FRF_UN_rip), 'r-', 'LineWidth', 2);
% hold on
% plot(omega, abs(FRF_UN), 'k--', 'LineWidth', 2);
% title('FRF senza Smorzamento');
% xlabel('Frequenza [Hz]');
% ylabel('|FRF|');
% grid on;
% 
% figure;
% plot(omega, abs(FRF_PROP_DUMP_rip), 'r-', 'LineWidth', 2);
% hold on
% plot(omega, abs(FRF_PROP_DUMP), 'k--', 'LineWidth', 2);
% title('FRF senza Smorzamento');
% xlabel('Frequenza [Hz]');
% ylabel('|FRF|');
% grid on;
% 
% figure;
% plot(omega, abs(FRF_NO_PROP_rip), 'r-', 'LineWidth', 2);
% hold on
% plot(omega, abs(FRF_NO_PROP), 'k--', 'LineWidth', 2);
% title('FRF senza Smorzamento');
% xlabel('Frequenza [Hz]');
% ylabel('|FRF|');
% grid on;