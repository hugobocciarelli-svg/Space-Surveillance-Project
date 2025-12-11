clear; clc; close all;

Target_selected = 21529 ; %id de la target

%% 1.Extraction des mesures et conversion en ECI
meas = Extract_data(Target_selected);
Z_ECEF=LatLonDist2ECEF(meas.azimuth,meas.elevation,meas.distance); %conversion des données az,elev,dist en ECEF
Z_eci=ECEFtoECI(Z_ECEF,meas.datetime);

idx_meas=1 ; %index de la mesure actuelle
Pos = [meas.x_true'; meas.y_true'; meas.z_true']; %vecteur des données réelles dans ECEF
%% 2. PARAMÈTRES DE SIMULATION
t0 = meas.time(1);                  % Temps initial [s]
dt = 1.6/2;                           % Pas de temps [s](Synchro sur le dt mesure pour des questions de simplicités)
t_sim = 3600*12;                      % Temps de simulation [s] 
t_end = t0 + t_sim;                 % Temps  final [s] 
TimeT = t0:dt:t_end;                % Vecteur temps
n_steps = length(TimeT);            % Nombre de pas

%% 3. CONDITIONS INITIALES
% État initial: [x, y, z, vx, vy, vz] en km et km/s
X0_real = [-4.287964709629751e+02; -4.788795050327217e+03; 5.344726876143979e+03; 0; 7.5; 0];                  % etat initial
X0_est = X0_real + [10; 5; -5; 0.1; -0.1; 0.05];    % Estimation initiale (avec erreur)
n_meas = size(Z_eci, 2);%Nombre de mesures prises en compte   
%% 4. PARAMÈTRES DU FILTRE DE KALMAN
% Matrice de mesure 
H = [1 0 0 0 0 0;
     0 1 0 0 0 0;
     0 0 1 0 0 0];

% Covariance initiale de l'estimation
P_est = diag([100, 100, 100, 0.1, 0.1, 0.1]); 

% Bruit du processus (incertitude du modèle)
W = diag([0.001, 0.001, 0.001, 0.00001, 0.00001, 0.00001]);

% Bruit de mesure (incertitude des capteurs)
%V = diag([0.1, 0.1, 0.1]); % Variance de 25 km² sur chaque position
V = diag([0, 0, 0]);
%% 5. INITIALISATION DES VARIABLES DE STOCKAGE
% Pour les vecteurs corrigés (+ initial)
Corrected_states = zeros(6, n_steps);
Corrected_states(:, 1) = X0_est;

% Pour les prédictions (avant correction)
predicted_states = zeros(6, n_steps);
predicted_states(:, 1) = X0_est;

% Pour la covariance
P_history = cell(1, n_steps);
P_history{1} = P_est;  
compteur_corr=0;

%% 6. BOUCLE PRINCIPALE - SIMULATION ET FILTRAGE
for i = 1:n_steps-1
    t = TimeT(i);
    
    %Prédiction du vecteur d'etat 
    x_pred = RungeKutta4(@modele_dynamique, t, Corrected_states(:, i), dt);
    predicted_states(:, i+1) = x_pred;

    % Matrice de transition d'état (Jacobienne)
    F=MatriceTransition(x_pred,dt);

    %Propagation de la covariance
    P_pred = F * P_history{i} * F' + W;
    
    % Cherche si une mesure existe proche du temps actuel (tolérance de dt/2)
    measure_dispo = false;
    if idx_meas <= n_meas
        time_diff = abs(meas.time(idx_meas) - t);
        if time_diff < dt/2  % Tolérance
            measure_dispo = true;
        end
    end
    
    if measure_dispo
        compteur_corr=compteur_corr+1;
        % extraction de la mesure actuelle 
        z_k = Z_eci(:, idx_meas);
        % Mise à jour de Kalman
        [X, P] = KalmanFilter(x_pred, P_pred, z_k, H, V);
        % Stocker les résultats corrigés
        Corrected_states(:, i+1) = X;
        P_history{i+1} = P;
        
        % Passer à la mesure suivante
        idx_meas = idx_meas + 1;
        
    else     %Si pas de mesure disponible on fait une estimation sans correction
        Corrected_states(:, i+1) = x_pred;
        P_history{i+1} = P_pred;
    end
end
delta = predicted_states-Corrected_states;

%% 7. Affichage

PlotTrajectoryComparison(TimeT, Corrected_states, Pos,Z_eci, meas.time, meas.datetime);

VisualizeOrbitAerospace(Corrected_states,TimeT,meas,t_sim)