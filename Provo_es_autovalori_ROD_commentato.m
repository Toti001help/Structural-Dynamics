clc;
clear;
close;

%% Informazioni sulla trave
L = 10;        % Lunghezza della trave (m)
E = 210e9;     % Modulo di Young (Pa)
I = 8.33e-6;   % Momento d'inerzia della sezione (m^4)
rho = 7800;    % Densità del materiale (kg/m^3)
A = 0.01;      % Area della sezione trasversale (m^2)

%% Soluzione analitica
N_eigen=9;  % Numero di autovalori da calcolare
% Vettore degli indici degli autovalori
a=linspace(1,N_eigen,N_eigen);
% Calcolo delle pulsazioni analitiche
omega_n_vec=a*sqrt(E/rho)*(pi/(L));

%% Soluzione numerica al variare del numero di elementi
l_bound=2;  % Limite inferiore per il numero di elementi
up_bound=36; % Limite superiore per il numero di elementi
N_vec=(l_bound:1:up_bound); % Vettore con il numero di elementi considerati
N_free_vec = N_vec-1; % Numero di gradi di libertà liberi
n_it=length(N_vec); % Numero di iterazioni
err=nan(n_it,N_eigen); % Matrice per l'errore relativo
omega_mat=nan(n_it,N_eigen); % Matrice per gli autovalori numerici

for i=1:n_it
    N=N_vec(i); % Numero di elementi
    N_free = N_free_vec(i); % Numero di GDL liberi
    L_el=L/N;  % Lunghezza di un elemento
    
    % Matrice di rigidezza per un elemento asta 
    k_el = (E * A / L_el) * [1, -1; -1, 1];
    
    % Matrice di massa per un elemento asta (formulazione consistente)
    m_el = (rho * A * L_el / 6) * [2, 1; 1, 2];
    
    % Inizializzazione delle matrici globali
    K = zeros(N+1, N+1);  % Matrice di rigidezza globale
    M = zeros(N+1, N+1);  % Matrice di massa globale
    
    % Assemblaggio delle matrici globali
    for p = 1:N
        nodes = [p, p+1];  % Nodi dell'elemento
        K(nodes, nodes) = K(nodes, nodes) + k_el;
        M(nodes, nodes) = M(nodes, nodes) + m_el;
    end
    
    % Applicazione delle condizioni al contorno
    free_nodes = 2:N; % Nodi liberi (asta vincolata agli estremi)
    
    % Riduzione del sistema ai soli nodi liberi
    K = K(free_nodes, free_nodes);    % Matrice di rigidezza ridotta
    M = M(free_nodes, free_nodes);    % Matrice di massa ridotta
    
    % Risoluzione del problema agli autovalori generalizzato
    [Phi, lambda_vec] = eig(K, M);
    lambda_vec=sort(diag(lambda_vec)); % Ordinamento degli autovalori
    lambda_vec=lambda_vec';
    
    % Determinazione del numero massimo di autovalori da considerare
    minim=min(N_eigen,floor((N_free+2)/2));
    
    % Memorizzazione delle pulsazioni numeriche e calcolo dell'errore
    omega_mat(i,1:minim) = sqrt(lambda_vec(1:minim));
    err(i,1:minim) = 100*abs((sqrt(lambda_vec(1:minim))-omega_n_vec(1:minim)))./omega_n_vec(1:minim);
end

%% Grafica dei risultati
r=20; % Numero di punti da rappresentare

% Creazione dei vettori delle pulsazioni analitiche
omega_1_vec=ones(1,r)*omega_n_vec(1);
omega_2_vec=ones(1,r)*omega_n_vec(2);
omega_3_vec=ones(1,r)*omega_n_vec(3);
omega_4_vec=ones(1,r)*omega_n_vec(4);

% Grafico degli autovalori analitici e numerici
figure(1)
plot(N_free_vec(1:r),omega_1_vec,N_free_vec(1:r),omega_2_vec,...
    N_free_vec(1:r),omega_3_vec,N_free_vec(1:r),omega_4_vec,...
    "LineStyle","-","LineWidth",1);
hold on
grid on
plot(N_free_vec(1:r),omega_mat((1:r),1),N_free_vec(1:r),omega_mat((1:r),2),...
    N_free_vec(1:r),omega_mat((1:r),3),N_free_vec(1:r),omega_mat((1:r),4),...
    "LineStyle",":","Marker",".","MarkerSize",16);

title('Autovalori per asta, appoggiata-appoggiata');
ylabel('$\omega_{n}(rad/s)$','Interpreter','latex','FontSize',12);
xlabel('Numero di GDL','FontSize',12);
lgd = legend('$\omega_{1,an}$','$\omega_{2,an}$','$\omega_{3,an}$','$\omega_{4,an}$',...
    '$\omega_{1,num}$','$\omega_{2,num}$','$\omega_{3,num}$','$\omega_{4,num}$');
lgd.Interpreter = 'latex'; 
lgd.FontSize = 11;

% Grafico dell'errore percentuale
figure(2)
plot(N_free_vec,err(:,1),N_free_vec,err(:,2),N_free_vec,err(:,3),N_free_vec,err(:,4),...
    N_free_vec,err(:,5),N_free_vec,err(:,6),N_free_vec,err(:,7),N_free_vec,err(:,8), ...
    N_free_vec,err(:,9),"LineStyle",":","Marker",".","MarkerSize",12)
hold on
grid on

title('Autovalori per asta, appoggiata-appoggiata');
ylabel('Err %','Interpreter','latex','FontSize',12);
xlabel('Numero di GDL','FontSize',12);
lgd = legend('$\omega_{1}$','$\omega_{2}$','$\omega_{3}$','$\omega_{4}$',...
    '$\omega_{5}$','$\omega_{6}$','$\omega_{7}$','$\omega_{8}$','$\omega_{9}$');
lgd.Interpreter = 'latex'; 
lgd.FontSize = 11;
