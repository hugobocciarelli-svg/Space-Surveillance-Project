%% ==================== PROJET SURVEILLANCE SPATIALE =====================
%  Détermination d'orbite et détection de collision avec perturbations
%  Auteurs: [Votre équipe]
%  Date: [Date]
%  Description: Filtre de Kalman étendu pour l'estimation et la prédiction
%               d'orbite avec perturbations J2, traînée et radiation solaire
% ========================================================================

clear; clc; close all;
tic; % Début du chronomètre

fprintf('========================================\n');
fprintf('   DÉTERMINATION D\'ORBITE AVEC PERTURBATIONS\n');
fprintf('========================================\n\n');

%% ==================== 1. CONFIGURATION ================================
% Choix de la cible à suivre
Target_selected = 22;      % ID de la cible (à changer selon besoin)

% Paramètres de simulation
dt = 1.6;                  % Pas de temps d'intégration [s]
t_sim = 3600;              % Temps de simulation [s] (1 heure)
fprintf('Configuration: Target=%d, dt=%.2fs, t_sim=%.0fs (%.1f heures)\n', ...
        Target_selected, dt, t_sim, t_sim/3600);

%% ==================== 2. PARAMÈTRES DU SATELLITE ======================
% --- Paramètres fondamentaux ---
sat.mu = 398600.4418;      % Paramètre gravitationnel Terre [km^3/s^2]
sat.R_terre = 6378.1363;   % Rayon terrestre équatorial [km]
sat.J2 = 1.08263e-3;       % Coefficient d'aplatissement J2 (unité: sans)

% --- Paramètres de TRAÎNÉE atmosphérique ---
sat.Cd = 2.2;              % Coefficient de traînée (typique 2.0-2.3)
sat.A_drag_m2 = 1.0;       % Surface de référence pour la traînée [m^2]
sat.masse_kg = 1000;       % Masse du satellite [kg]
% Coefficient balistique (conversion d'unités: m²/kg → km²/kg)
sat.BC = (sat.Cd * sat.A_drag_m2 / sat.masse_kg) * 1e6;

% --- Paramètres de RADIATION SOLAIRE ---
sat.Cr = 1.3;              % Coefficient de réflexion (1: absorption, 2: réflexion)
sat.A_srp_m2 = 1.0;        % Surface exposée au soleil [m^2]
% Pression solaire à 1 UA [N/m² = kg/(m·s²)]
P_sol = 4.56e-6;
% Coefficient de radiation solaire [km²/kg]
sat.K_srp = P_sol * sat.Cr * (sat.A_srp_m2 / sat.masse_kg) * 1e6;

% --- Paramètres du modèle de densité atmosphérique exponentiel ---
% Modèle: rho = rho0 * exp(-(altitude - alt0) / H)
sat.rho0 = 1.225e9;        % Densité au niveau de la mer [kg/km³] (équiv. à 1.225 kg/m³)
sat.H = 8.5;               % Hauteur d'échelle [km]

% --- Rotation terrestre pour la traînée ---
sat.omega_terre = 7.2921150e-5;  % Vitesse angulaire Terre [rad/s]

fprintf('Paramètres satellite chargés:\n');
fprintf('  - Masse: %.0f kg, Surface traînée: %.2f m²\n', sat.masse_kg, sat.A_drag_m2);
fprintf('  - BC: %.4f km²/kg, K_srp: %.2e km²/kg\n', sat.BC, sat.K_srp);

%% ==================== 3. EXTRACTION DES MESURES =======================
fprintf('\n--- Extraction des données pour la cible %d ---\n', Target_selected);

% Extraction avec VOTRE fonction
Target_data = Extract_data(Target_selected);

% Préparation des vecteurs de mesures en ECEF
Z_ecef = [Target_data.x_true'; Target_data.y_true'; Target_data.z_true'];
n_meas = size(Z_ecef, 2);

% Conversion des mesures ECEF → ECI
fprintf('Conversion ECEF → ECI des %d mesures...\n', n_meas);
Z_eci = zeros(3, n_meas);

for k = 1:n_meas
    r_ecef_k = Z_ecef(:, k);
    utc_datetime_k = Target_data.datetime(k);
    
    % Conversion ECEF → ECI (votre fonction)
    [r_eci_k, ~] = ecef2eci(utc_datetime_k, r_ecef_k);
    Z_eci(:, k) = r_eci_k;
    
    % Barre de progression
    if mod(k, floor(n_meas/10)) == 0
        fprintf('  %.0f%%\n', (k/n_meas)*100);
    end
end
fprintf('Conversion terminée.\n');

idx_meas = 1;  % Index pour parcourir les mesures

%% ==================== 4. PARAMÈTRES TEMPORELS =========================
t0 = Target_data.time(1);       % Temps initial [s] (première mesure)
t_end = t0 + t_sim;             % Temps final de simulation [s]
TimeT = t0:dt:t_end;            % Vecteur temps
n_steps = length(TimeT);        % Nombre de pas de simulation

fprintf('\nSimulation temporelle:\n');
fprintf('  t0: %.1f s (première mesure)\n', t0);
fprintf('  t_end: %.1f s\n', t_end);
fprintf('  Nombre de pas: %d (dt=%.2fs)\n', n_steps, dt);

%% ==================== 5. CONDITIONS INITIALES =========================
% État réel initial (inconnu en pratique - pour validation seulement)
% Orbite circulaire à 7078 km d'altitude
X0_real = [7078; 0; 0; 0; 7.5; 0];  % [x,y,z,vx,vy,vz] en km et km/s

% Estimation initiale (avec erreur délibérée)
X0_est = X0_real + [10; 5; -5; 0.1; -0.1; 0.05];

fprintf('\nConditions initiales:\n');
fprintf('  Vérité: [%.1f, %.1f, %.1f, %.2f, %.2f, %.2f]\n', X0_real);
fprintf('  Estimation: [%.1f, %.1f, %.1f, %.2f, %.2f, %.2f]\n', X0_est);

%% ==================== 6. PARAMÈTRES DU FILTRE DE KALMAN ===============
fprintf('\n--- Configuration du filtre de Kalman étendu ---\n');

% Matrice de mesure (on mesure uniquement la position)
H = [1 0 0 0 0 0;
     0 1 0 0 0 0;
     0 0 1 0 0 0];

% Covariance initiale de l'estimation
% Réflète notre confiance dans l'estimation initiale
P_est = diag([100, 100, 100, 0.1, 0.1, 0.1]);  % [km², km², km², (km/s)², ...]
fprintf('  Covariance initiale diagonale:\n');
fprintf('    Position: σ=%.1f km\n', sqrt(P_est(1,1)));
fprintf('    Vitesse: σ=%.3f km/s\n', sqrt(P_est(4,4)));

% Bruit du processus (incertitude du modèle dynamique)
% Représente les erreurs de modélisation des perturbations
W = diag([0.001, 0.001, 0.001, 0.00001, 0.00001, 0.00001]);
fprintf('  Bruit processus: Q=diag(%s)\n', mat2str(diag(W)', 3));

% Bruit de mesure (incertitude des capteurs)
% Correspond à la précision des mesures de position
V = diag([2, 2, 2]);  % Variance de 2 km² sur chaque position
fprintf('  Bruit mesure: R=diag(%s) (σ=%.2f km)\n', mat2str(diag(V)', 2), sqrt(V(1,1)));

%% ==================== 7. INITIALISATION ===============================
fprintf('\n--- Initialisation des variables de stockage ---\n');

% Pour les états corrigés (après mise à jour Kalman)
Corrected_states = zeros(6, n_steps);
Corrected_states(:, 1) = X0_est;

% Pour les états prédits (avant correction)
predicted_states = zeros(6, n_steps);
predicted_states(:, 1) = X0_est;

% Pour les covariances
P_history = cell(1, n_steps);
P_history{1} = P_est;

% Pour la vérité terrain (validation)
Truth_states = zeros(6, n_steps);
Truth_states(:, 1) = X0_real;

% Compteurs
compteur_corr = 0;        % Nombre de corrections appliquées
compteur_pred = 0;        % Nombre de prédictions sans mesure

% Préparation du modèle dynamique avec paramètres
fprintf('  Préparation du modèle dynamique avec perturbations...\n');
modele_dyn_complet = @(t, Y) modele_dynamique(t, Y, sat);

%% ==================== 8. BOUCLE PRINCIPALE ============================
fprintf('\n========================================\n');
fprintf('   DÉBUT DE LA SIMULATION - FILTRAGE\n');
fprintf('========================================\n\n');

for i = 1:n_steps-1
    t = TimeT(i);
    
    % Affichage de progression
    if mod(i, floor(n_steps/10)) == 0
        fprintf('  Progression: %.0f%% (t=%.1f s)\n', (i/n_steps)*100, t);
    end
    
    % === ÉTAPE DE PRÉDICTION (propagation de l'état) ===
    x_pred = RungeKutta4(modele_dyn_complet, t, Corrected_states(:, i), dt);
    predicted_states(:, i+1) = x_pred;
    compteur_pred = compteur_pred + 1;
    
    % === ÉTAPE DE PRÉDICTION (propagation de la covariance) ===
    % Calcul de la matrice de transition d'état (STM)
    F = MatriceTransition(x_pred, dt, sat);
    
    % Discrétisation du bruit de processus
    % Approximation: Q = W * dt (pour des petits dt)
    Q = W * dt;
    
    % Propagation de la covariance: P_pred = F * P * F' + Q
    P_pred = F * P_history{i} * F' + Q;
    
    % === VÉRIFICATION DISPONIBILITÉ MESURE ===
    measure_dispo = false;
    if idx_meas <= n_meas
        % On cherche une mesure proche du temps courant (tolérance dt/2)
        time_diff = abs(Target_data.time(idx_meas) - t);
        if time_diff < dt/2
            measure_dispo = true;
        end
    end
    
    % === ÉTAPE DE CORRECTION (si mesure disponible) ===
    if measure_dispo
        compteur_corr = compteur_corr + 1;
        
        % Extraction de la mesure courante
        z_k = Z_eci(:, idx_meas);
        
        % Application du filtre de Kalman
        [x_updated, P_updated] = KalmanFilter(x_pred, P_pred, z_k, H, V);
        
        % Stockage des résultats
        Corrected_states(:, i+1) = x_updated;
        P_history{i+1} = P_updated;
        
        % Passage à la mesure suivante
        idx_meas = idx_meas + 1;
        
    else
        % Pas de mesure disponible → l'état prédit devient l'état corrigé
        Corrected_states(:, i+1) = x_pred;
        P_history{i+1} = P_pred;
    end
    
    % === PROPAGATION DE LA VÉRITÉ TERRAIN (pour validation) ===
    Truth_states(:, i+1) = RungeKutta4(modele_dyn_complet, t, Truth_states(:, i), dt);
end

fprintf('\n--- Simulation terminée ---\n');
fprintf('  Prédictions: %d\n', compteur_pred);
fprintf('  Corrections: %d (%.1f%%)\n', compteur_corr, (compteur_corr/compteur_pred)*100);

%% ==================== 9. ANALYSE DES RÉSULTATS ========================
fprintf('\n========================================\n');
fprintf('   ANALYSE DES RÉSULTATS\n');
fprintf('========================================\n\n');

% Calcul des erreurs d'estimation
pos_err = vecnorm(Corrected_states(1:3,:) - Truth_states(1:3,:), 2, 1);
vel_err = vecnorm(Corrected_states(4:6,:) - Truth_states(4:6,:), 2, 1);

% Statistiques des erreurs finales
fprintf('Erreurs finales (dernier pas):\n');
fprintf('  Position: %.3f km (cible: < 1.0 km)\n', pos_err(end));
fprintf('  Vitesse: %.5f km/s (cible: < 0.01 km/s)\n', vel_err(end));

% Évolution des incertitudes (écarts-types)
pos_std = zeros(3, n_steps);
for i = 1:n_steps
    P = P_history{i};
    pos_std(:, i) = sqrt(diag(P(1:3, 1:3)));
end

fprintf('\nIncertitudes moyennes (1σ):\n');
fprintf('  σ_x: %.2f km, σ_y: %.2f km, σ_z: %.2f km\n', ...
        mean(pos_std(1,:)), mean(pos_std(2,:)), mean(pos_std(3,:)));

%% ==================== 10. VISUALISATIONS ==============================
fprintf('\n--- Génération des graphiques ---\n');

% Figure 1: Évolution des erreurs
figure('Position', [100, 100, 1200, 500], 'Name', 'Erreurs d''estimation');

subplot(1,2,1);
plot(TimeT, pos_err, 'b-', 'LineWidth', 1.5);
xlabel('Temps [s]', 'FontSize', 11);
ylabel('Erreur de position [km]', 'FontSize', 11);
title('Erreur de position vs Temps', 'FontSize', 12);
grid on; box on;
xlim([TimeT(1), TimeT(end)]);
% Ligne horizontale pour la précision cible
hold on; plot([TimeT(1), TimeT(end)], [1.0, 1.0], 'r--', 'LineWidth', 1);
legend('Erreur réelle', 'Cible (1 km)', 'Location', 'best');

subplot(1,2,2);
plot(TimeT, vel_err, 'r-', 'LineWidth', 1.5);
xlabel('Temps [s]', 'FontSize', 11);
ylabel('Erreur de vitesse [km/s]', 'FontSize', 11);
title('Erreur de vitesse vs Temps', 'FontSize', 12);
grid on; box on;
xlim([TimeT(1), TimeT(end)]);
% Ligne horizontale pour la précision cible
hold on; plot([TimeT(1), TimeT(end)], [0.01, 0.01], 'r--', 'LineWidth', 1);
legend('Erreur réelle', 'Cible (0.01 km/s)', 'Location', 'best');

% Figure 2: Évolution des incertitudes
figure('Position', [100, 600, 1200, 400], 'Name', 'Incertitudes (écarts-types)');

subplot(1,3,1);
plot(TimeT, pos_std(1,:), 'b-', 'LineWidth', 1.5);
xlabel('Temps [s]'); ylabel('σ_x [km]');
title('Incertitude sur X'); grid on; box on;

subplot(1,3,2);
plot(TimeT, pos_std(2,:), 'g-', 'LineWidth', 1.5);
xlabel('Temps [s]'); ylabel('σ_y [km]');
title('Incertitude sur Y'); grid on; box on;

subplot(1,3,3);
plot(TimeT, pos_std(3,:), 'r-', 'LineWidth', 1.5);
xlabel('Temps [s]'); ylabel('σ_z [km]');
title('Incertitude sur Z'); grid on; box on;

% Figure 3: Trajectoire 3D
figure('Position', [700, 100, 800, 800], 'Name', 'Trajectoires 3D');

% Terre
[Xe, Ye, Ze] = sphere(50);
R_T = sat.R_terre;
surf(Xe*R_T, Ye*R_T, Ze*R_T, ...
     'FaceColor', [0.2, 0.4, 0.9], ...
     'EdgeColor', 'none', ...
     'FaceAlpha', 0.7);
hold on; grid on; axis equal;

% Trajectoire estimée
plot3(Corrected_states(1,:), Corrected_states(2,:), Corrected_states(3,:), ...
      'b-', 'LineWidth', 2, 'DisplayName', 'Orbite estimée');

% Trajectoire réelle
plot3(Truth_states(1,:), Truth_states(2,:), Truth_states(3,:), ...
      'g--', 'LineWidth', 1.5, 'DisplayName', 'Orbite réelle');

% Points de mesures (convertis en ECI pour l'affichage)
plot3(Z_eci(1,:), Z_eci(2,:), Z_eci(3,:), ...
      'ro', 'MarkerSize', 4, 'MarkerFaceColor', 'r', ...
      'DisplayName', 'Mesures');

% Point de départ
plot3(X0_est(1), X0_est(2), X0_est(3), ...
      'k^', 'MarkerSize', 10, 'MarkerFaceColor', 'y', ...
      'DisplayName', 'Départ (estimé)');

xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('Trajectoire orbitale', 'FontSize', 14);
legend('Location', 'best'); view(45, 25);

% Figure 4: Énergie orbitale (conservation)
figure('Position', [1300, 100, 800, 400], 'Name', 'Énergie orbitale');

% Calcul de l'énergie spécifique
energy_est = zeros(1, n_steps);
energy_true = zeros(1, n_steps);

for i = 1:n_steps
    r_est = norm(Corrected_states(1:3,i));
    v_est = norm(Corrected_states(4:6,i));
    energy_est(i) = v_est^2/2 - sat.mu/r_est;
    
    r_true = norm(Truth_states(1:3,i));
    v_true = norm(Truth_states(4:6,i));
    energy_true(i) = v_true^2/2 - sat.mu/r_true;
end

plot(TimeT, energy_est, 'b-', 'LineWidth', 1.5); hold on;
plot(TimeT, energy_true, 'g--', 'LineWidth', 1.5);
xlabel('Temps [s]'); ylabel('Énergie spécifique [km²/s²]');
title('Énergie orbitale (devrait être ~constante)');
legend('Estimée', 'Réelle'); grid on; box on;

%% ==================== 11. RAPPORT DE PERFORMANCE ======================
fprintf('\n========================================\n');
fprintf('   RAPPORT DE PERFORMANCE\n');
fprintf('========================================\n\n');

% Calcul des moyennes et écarts-types
mean_pos_err = mean(pos_err);
std_pos_err = std(pos_err);
max_pos_err = max(pos_err);

mean_vel_err = mean(vel_err);
std_vel_err = std(vel_err);
max_vel_err = max(vel_err);

fprintf('ERREUR DE POSITION:\n');
fprintf('  Moyenne: %.3f km\n', mean_pos_err);
fprintf('  Écart-type: %.3f km\n', std_pos_err);
fprintf('  Maximum: %.3f km\n', max_pos_err);
fprintf('  → Précision relative: %.2f%%\n', (mean_pos_err/7078)*100);

fprintf('\nERREUR DE VITESSE:\n');
fprintf('  Moyenne: %.5f km/s\n', mean_vel_err);
fprintf('  Écart-type: %.5f km/s\n', std_vel_err);
fprintf('  Maximum: %.5f km/s\n', max_vel_err);
fprintf('  → Précision relative: %.2f%%\n', (mean_vel_err/7.5)*100);

fprintf('\nCONSERVATION DE L''ÉNERGIE:\n');
fprintf('  Variation énergie estimée: %.4f km²/s²\n', max(energy_est)-min(energy_est));
fprintf('  Variation énergie réelle: %.4f km²/s²\n', max(energy_true)-min(energy_true));

% Évaluation binaire des objectifs
fprintf('\nOBJECTIFS ATTEINTS:\n');
if mean_pos_err < 1.0
    fprintf('  ✓ Précision position < 1.0 km\n');
else
    fprintf('  ✗ Précision position > 1.0 km (%.2f km)\n', mean_pos_err);
end

if mean_vel_err < 0.01
    fprintf('  ✓ Précision vitesse < 0.01 km/s\n');
else
    fprintf('  ✗ Précision vitesse > 0.01 km/s (%.3f km/s)\n', mean_vel_err);
end

% Temps d'exécution
exec_time = toc;
fprintf('\nTEMPS D''EXÉCUTION: %.2f secondes\n', exec_time);

%% ==================== 12. SAUVEGARDE DES RÉSULTATS ====================
fprintf('\n--- Sauvegarde des résultats ---\n');

% Structure de résultats
results = struct();
results.Time = TimeT;
results.Truth = Truth_states;
results.Estimated = Corrected_states;
results.Predicted = predicted_states;
results.P_history = P_history;
results.Errors.pos = pos_err;
results.Errors.vel = vel_err;
results.Errors.energy_est = energy_est;
results.Errors.energy_true = energy_true;
results.SatParams = sat;
results.KalmanParams.H = H;
results.KalmanParams.W = W;
results.KalmanParams.V = V;
results.Simulation.dt = dt;
results.Simulation.n_steps = n_steps;
results.Simulation.n_corrections = compteur_corr;

% Sauvegarde dans un fichier .mat
save('results_orbit_determination.mat', 'results');
fprintf('  Résultats sauvegardés dans: results_orbit_determination.mat\n');

% Export CSV pour analyse externe
csv_data = [TimeT', Corrected_states', Truth_states', pos_err', vel_err'];
csv_header = {'Time', 'X_est', 'Y_est', 'Z_est', 'Vx_est', 'Vy_est', 'Vz_est', ...
              'X_true', 'Y_true', 'Z_true', 'Vx_true', 'Vy_true', 'Vz_true', ...
              'Pos_Err', 'Vel_Err'};
fid = fopen('orbit_results.csv', 'w');
fprintf(fid, '%s,', csv_header{1:end-1});
fprintf(fid, '%s\n', csv_header{end});
fclose(fid);
dlmwrite('orbit_results.csv', csv_data, '-append', 'precision', '%.6f');
fprintf('  Données exportées dans: orbit_results.csv\n');

%% ==================== 13. PRÉPARATION POUR LA DÉTECTION ===============
fprintf('\n========================================\n');
fprintf('   PRÉPARATION POUR LA DÉTECTION DE COLLISION\n');
fprintf('========================================\n\n');

% Extraction de l'état final et de sa covariance
X_final = Corrected_states(:, end);
P_final = P_history{end};

fprintf('État final estimé:\n');
fprintf('  Position: [%.1f, %.1f, %.1f] km\n', X_final(1:3));
fprintf('  Vitesse: [%.3f, %.3f, %.3f] km/s\n', X_final(4:6));

fprintf('\nPour la détection de collision, vous aurez besoin de:\n');
fprintf('  1. Propager cet état dans le futur (fonction RK4)\n');
fprintf('  2. Pour chaque satellite, avoir X(t) et P(t)\n');
fprintf('  3. Calculer la distance relative: Δr = r_A - r_B\n');
fprintf('  4. Calculer la covariance relative: P_rel = P_A + P_B\n');
fprintf('  5. Appliquer le critère de Patera (2001) pour la probabilité\n');

fprintf('\n========================================\n');
fprintf('   SIMULATION TERMINÉE AVEC SUCCÈS\n');
fprintf('========================================\n');

%% ==================== FONCTIONS LOCALES ===============================
% Note: Les fonctions suivantes doivent être dans des fichiers séparés .m
% Elles sont décrites ici pour référence

% function dot_etat = modele_dynamique(t, etat, sat)
% function Y_next = RungeKutta4(model, t, Y, dt)
% function F = MatriceTransition(x, dt, sat)
% function [x_updated, P_updated] = KalmanFilter(x_pred, P_pred, Z, H, V)
% function [r_eci, v_eci] = ecef2eci(utc_datetime, r_ecef, v_ecef)
% function [Target_data] = Extract_data(n)

%% ==================== FIN DU SCRIPT ===================================