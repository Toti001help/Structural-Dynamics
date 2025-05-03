clc; 
clear;
close all;

format long  % Per visualizzare più cifre decimali nel workspace

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

for i = 1:size(PHI_UN, 2)
    v = PHI_UN(:, i);
    norm_factor = sqrt(v' * M * v);
    PHI_UN(:, i) = v / norm_factor;
end

% Forma polare degli autovettori
PHI_UN_POL = strings(size(PHI_UN));
for i = 1:size(PHI_UN,1)
    for j = 1:size(PHI_UN,2)
        modulo = abs(PHI_UN(i,j));
        fase = rad2deg(angle(PHI_UN(i,j)));
        PHI_UN_POL(i,j) = sprintf('%.6f ∠ %.4f°', modulo, fase);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%        PROPORTIONAL STRUCTURAL DAMPING               %%%%%%%%%%%%%

beta=0.05;
D_1 = beta*K;
[PHI_PROP_DUMP, LAMBDA_PROP_DUMP] = eig(K+1j*D_1, M);
[LAMBDA_PROP_DUMP,ind_2]=sort(diag(LAMBDA_PROP_DUMP));
LAMBDA_PROP_DUMP=diag(LAMBDA_PROP_DUMP);
PHI_PROP_DUMP=PHI_PROP_DUMP(:,ind_2);

for i = 1:size(PHI_PROP_DUMP, 2)
    v = PHI_PROP_DUMP(:, i);
    norm_factor = sqrt(conj(v)' * M * v);  % Normalizzazione Hermitiana
    PHI_PROP_DUMP(:, i) = v / norm_factor;
end

% Forma polare degli autovettori
PHI_PROP_POL = strings(size(PHI_PROP_DUMP));
for i = 1:size(PHI_PROP_DUMP,1)
    for j = 1:size(PHI_PROP_DUMP,2)
        modulo = abs(PHI_PROP_DUMP(i,j));
        fase = rad2deg(angle(PHI_PROP_DUMP(i,j)));
        PHI_PROP_POL(i,j) = sprintf('%.3f ∠ %.1f°', modulo, fase);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%        NON PROPORTIONAL STRUCTURAL DAMPING           %%%%%%%%%%%%%

d=0.3*k;
D_2 = zeros(3,3);
D_2(1,1) = d;
[PHI_NOPROP, LAMBDA_NO_PROP] = eig(K+(1j*D_2), M);
[LAMBDA_NO_PROP,ind_3]=sort(diag(LAMBDA_NO_PROP));
LAMBDA_NO_PROP=diag(LAMBDA_NO_PROP);
PHI_NOPROP=PHI_NOPROP(:,ind_3);

for i = 1:size(PHI_NOPROP, 2)
    v = PHI_NOPROP(:, i);
    norm_factor = sqrt(conj(v)' * M * v);
    PHI_NOPROP(:, i) = v / norm_factor;
end

% Forma polare degli autovettori
PHI_NOPROP_POL = strings(size(PHI_NOPROP));
for i = 1:size(PHI_NOPROP,1)
    for j = 1:size(PHI_NOPROP,2)
        modulo = abs(PHI_NOPROP(i,j));
        fase = rad2deg(angle(PHI_NOPROP(i,j)));
        PHI_NOPROP_POL(i,j) = sprintf('%.3f ∠ %.1f°', modulo, fase);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%     DISPLAY (chiaro in forma matriciale)     %%%%%%%%%%

fprintf('\n--- AUTOVETTORI NON SMORZATI (forma polare) ---\n');
for i = 1:size(PHI_UN_POL, 1)
    fprintf('%25s %25s %25s\n', PHI_UN_POL(i,1), PHI_UN_POL(i,2), PHI_UN_POL(i,3));
end

fprintf('\n--- AUTOVETTORI CON SMORZAMENTO PROPORZIONALE (forma polare) ---\n');
for i = 1:size(PHI_PROP_POL, 1)
    fprintf('%25s %25s %25s\n', PHI_PROP_POL(i,1), PHI_PROP_POL(i,2), PHI_PROP_POL(i,3));
end

fprintf('\n--- AUTOVETTORI CON SMORZAMENTO NON PROPORZIONALE (forma polare) ---\n');
for i = 1:size(PHI_NOPROP_POL, 1)
    fprintf('%25s %25s %25s\n', PHI_NOPROP_POL(i,1), PHI_NOPROP_POL(i,2), PHI_NOPROP_POL(i,3));
end


% grafica


% Supponiamo che tu abbia già calcolato gli autovettori PHI_UN, PHI_PROP_DUMP, PHI_NOPROP

% Numero di modi (m) e DOF
modes = 3;  % Assumiamo che ci siano 3 modi
DOF = 3;  % 3 gradi di libertà

% Grafico degli autovettori nel piano complesso
figure;

% Caso senza smorzamento
subplot(3, 1, 1);
hold on;
for i = 1:modes
    % Estrai l'autovettore i-esimo per il caso senza smorzamento
    v = PHI_UN(:, i);
    
    % Traccia nel piano complesso
    plot(real(v), imag(v), 'o-', 'LineWidth', 2, 'MarkerSize', 8);
end
title('Autovettori senza Smorzamento');
xlabel('Parte Reale');
ylabel('Parte Immaginaria');
legend('Modo 1', 'Modo 2', 'Modo 3');
grid on;

% Caso con smorzamento proporzionale
subplot(3, 1, 2);
hold on;
for i = 1:modes
    % Estrai l'autovettore i-esimo per il caso con smorzamento proporzionale
    v = PHI_PROP_DUMP(:, i);
    
    % Traccia nel piano complesso
    plot(real(v), imag(v), 'o-', 'LineWidth', 2, 'MarkerSize', 8);
end
title('Autovettori con Smorzamento Proporzionale');
xlabel('Parte Reale');
ylabel('Parte Immaginaria');
legend('Modo 1', 'Modo 2', 'Modo 3');
grid on;

% Caso con smorzamento non proporzionale
subplot(3, 1, 3);
hold on;
for i = 1:modes
    % Estrai l'autovettore i-esimo per il caso con smorzamento non proporzionale
    v = PHI_NOPROP(:, i);
    
    % Traccia nel piano complesso
    plot(real(v), imag(v), 'o-', 'LineWidth', 2, 'MarkerSize', 8);
end
title('Autovettori con Smorzamento Non Proporzionale');
xlabel('Parte Reale');
ylabel('Parte Immaginaria');
legend('Modo 1', 'Modo 2', 'Modo 3');
grid on;

% Miglioramenti grafici
sgtitle('Autovettori nel Piano Complesso per i Vari Tipi di Smorzamento', 'FontSize', 14);
