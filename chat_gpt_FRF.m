clc;
clear;
close all;

% Definizione delle matrici di massa, rigidezza e smorzamento (già note)
m1=1; % Kg
m2=0.95; % Kg
m3=1.05; % Kg
k=1.0e3; % N/m

M = [m1, 0, 0;
     0, m2, 0;
     0, 0, m3];

K = [k+k+k, -k, -k;
    -k, k+k+k, -k;
    -k, -k, k+k+k];

% Calcolare la frequenza naturale
[PHI_UN, LAMBDA_UN] = eig(K, M);
LAMBDA_UN = diag(LAMBDA_UN);
frequenze_naturali = sqrt(LAMBDA_UN);

% Definizione della frequenza
f1 = frequenze_naturali(1); % Frequenza naturale del primo modo
f3 = frequenze_naturali(3); % Frequenza naturale del terzo modo
omega1 = 2 * pi * f1;  % Frequenza angolare
omega3 = 2 * pi * f3;  % Frequenza angolare

% Range di frequenze per il calcolo della FRF
frequenze = linspace(1, 100, 3001);
omega = 2 * pi * frequenze;

% Forza applicata nel nodo 1 (in effetti una forza unitaria)
F = [1; 0; 0]; % Applicata sul primo nodo

% Risposta nel nodo 3
risposta = zeros(length(omega), 1);

% Calcolo della FRF per i tre modelli: senza smorzamento, proporzionale, non proporzionale

% Caso senza smorzamento
D_1 = zeros(3, 3);  % Matrice di smorzamento nulla
for i = 1:length(omega)
    % Matrice di trasferimento per la frequenza omega(i)
    A = M * (omega(i)^2 * eye(3) - K) + 1j * omega(i) * D_1;
    risposta(i) = F' * (A \ F);  % Calcola la risposta
end
FRF_UN = risposta;  % FRF senza smorzamento

% Caso con smorzamento proporzionale
beta = 0.05;  % Fattore di smorzamento
D_2 = beta * K;  % Matrice di smorzamento proporzionale
for i = 1:length(omega)
    % Matrice di trasferimento per la frequenza omega(i)
    A = M * (omega(i)^2 * eye(3) - K) + 1j * omega(i) * D_2;
    risposta(i) = F' * (A \ F);  % Calcola la risposta
end
FRF_PROP_DUMP = risposta;  % FRF con smorzamento proporzionale

% Caso con smorzamento non proporzionale
d = 0.3 * k;  % Smorzamento per il caso non proporzionale
D_3 = zeros(3, 3);
D_3(1, 1) = d;  % Assegnazione dello smorzamento non proporzionale
for i = 1:length(omega)
    % Matrice di trasferimento per la frequenza omega(i)
    A = M * (omega(i)^2 * eye(3) - K) + 1j * omega(i) * D_3;
    risposta(i) = F' * (A \ F);  % Calcola la risposta
end
FRF_NOPROP = risposta;  % FRF con smorzamento non proporzionale

% Plot delle FRF per i tre modelli
figure;

% FRF senza smorzamento
subplot(3, 1, 1);
plot(frequenze, abs(FRF_UN), 'r-', 'LineWidth', 2);
title('FRF senza Smorzamento');
xlabel('Frequenza [Hz]');
ylabel('|FRF|');
grid on;

% FRF con smorzamento proporzionale
subplot(3, 1, 2);
plot(frequenze, abs(FRF_PROP_DUMP), 'g-', 'LineWidth', 2);
title('FRF con Smorzamento Proporzionale');
xlabel('Frequenza [Hz]');
ylabel('|FRF|');
grid on;

% FRF con smorzamento non proporzionale
subplot(3, 1, 3);
plot(frequenze, abs(FRF_NOPROP), 'b-', 'LineWidth', 2);
title('FRF con Smorzamento Non Proporzionale');
xlabel('Frequenza [Hz]');
ylabel('|FRF|');
grid on;

% Miglioramenti grafici
sgtitle('Funzione di Risposta in Frequenza (FRF) per i Vari Modelli di Smorzamento', 'FontSize', 14);
