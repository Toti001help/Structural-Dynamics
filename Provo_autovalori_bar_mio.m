clc;
clear;
close;

%% informazioni sulla trave
L = 10;        % Lunghezza della trave (m)
E = 210e9;     % Modulo di Young (Pa)
I = 8.33e-6;   % Momento d'inerzia della sezione (m^4)
rho = 7800;    % Densità del materiale (kg/m^3)
A = 0.01;      % Area della sezione trasversale (m^2)

%% sol analitica
N_eigen=9;
a=linspace(1,N_eigen,N_eigen);
%sol trave appoggiata appoggiata
omega_n_vec=((a*pi/L).^2)*sqrt(E*I/(rho*A));

%% sol numerica al variare del numero di elementi
l_bound=1;
up_bound=17;
N_vec=(l_bound:1:up_bound);
N_free_vec = 2*N_vec;
n_it=length(N_vec);
err=nan(n_it,N_eigen); % ogni colonna errore relativo a fissato a utovalore
omega_mat=nan(n_it,N_eigen); % ogni colonna autovalore numerico relativo ...
                                % all'autovalore analitico



 for i=1:n_it
N_free = N_free_vec(i);
N=N_vec(i); %% N elementi, N+1 nodi

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

%lambda_vec=eig(M\K);
 [Phi, lambda_vec] = eig(K, M);
 lambda_vec=sort(diag(lambda_vec));
 lambda_vec=lambda_vec';
minim=min(N_eigen,floor((N_free+2)/2));
omega_mat(i,1:minim) = sqrt(lambda_vec(1:minim));
err(i,1:minim) = 100*abs((sqrt(lambda_vec(1:minim))-omega_n_vec(1:minim)))./omega_n_vec(1:minim);
end




%% grafica
% autovalori numerici e analitici
N_g=10;
omega_1_vec=ones(1,N_g)*omega_n_vec(1);
omega_2_vec=ones(1,N_g)*omega_n_vec(2);
omega_3_vec=ones(1,N_g)*omega_n_vec(3);
omega_4_vec=ones(1,N_g)*omega_n_vec(4);

figure(1)
plot(N_free_vec(1:N_g),omega_1_vec,N_free_vec(1:N_g),omega_2_vec,...
    N_free_vec(1:N_g),omega_3_vec,N_free_vec(1:N_g),omega_4_vec,...
    "LineStyle","-","LineWidth",1);%
hold on
grid on
plot(N_free_vec(1:N_g),omega_mat((1:N_g),1),N_free_vec(1:N_g),omega_mat((1:N_g),2),...
    N_free_vec(1:N_g),omega_mat((1:N_g),3),N_free_vec(1:N_g),omega_mat((1:N_g),4),...
    "LineStyle",":","Marker",".","MarkerSize",16);

    


    title('Autovalori per trave,appoggiata-appoggiata');
ylabel('$\omega_{n}(rad/s)$','Interpreter','latex','FontSize',12);
xlabel('Numero di GDL','FontSize',12);
lgd = legend('$\omega_{1,an}$','$\omega_{2,an}$','$\omega_{3,an}$','$\omega_{4,an}$',...
    '$\omega_{1,num}$','$\omega_{2,num}$','$\omega_{3,num}$','$\omega_{4,num}$');
lgd.Interpreter = 'latex'; 
lgd.FontSize = 11;

% errore percentuale
figure(2)
plot(N_free_vec,err(:,1),N_free_vec,err(:,2),N_free_vec,err(:,3),N_free_vec,err(:,4),...
    N_free_vec,err(:,5),N_free_vec,err(:,6),N_free_vec,err(:,7),N_free_vec,err(:,8), ...
    N_free_vec,err(:,9),"LineStyle",":","Marker",".","MarkerSize",12)
hold on
grid on

    title('Autovalori per trave, appoggiata-appoggiata');
ylabel('Err %','Interpreter','latex','FontSize',12);
xlabel('Numero di GDL','FontSize',12);
lgd = legend('$\omega_{1}$','$\omega_{2}$','$\omega_{3}$','$\omega_{4}$',...
    '$\omega_{5}$','$\omega_{6}$','$\omega_{7}$','$\omega_{8}$','$\omega_{9}$');
lgd.Interpreter = 'latex'; 
lgd.FontSize = 11;