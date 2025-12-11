% Script de traitement de données et formatage des données d'entrées
% Auteur: Hugo Bocciarelli
function [Target_data] = Extract_data(n)
    % Chargement des données
    loaded_data = load('mesures.mat');
    mesures = loaded_data.mesures;
    
    % Extraction des colonnes
    target_id = mesures(:, 1);
    time = mesures(:, 2);
    x_true = mesures(:, 3)/1000;
    y_true = mesures(:, 4)/1000;
    z_true = mesures(:, 5)/1000;
    azimuth = mesures(:, 6);
    elevation = mesures(:, 7);
    distance = mesures(:, 8)/1000;
    
    % Conversion du temps
    epoch_2000 = datetime(2000, 1, 1, 0, 0, 0);
    datetime_vals = epoch_2000 + seconds(time);
    
    % Affichage résumé
    fprintf('\n=== RÉSUMÉ DES DONNÉES ===\n');
    fprintf('Nombre de mesures: %d\n', size(mesures, 1));
    fprintf('Nombre de cibles: %d\n', length(unique(target_id)));
    fprintf('Période: %s à %s\n', datestr(datetime_vals(1)), datestr(datetime_vals(end)));
    fprintf('Plage azimuth: [%.2f, %.2f] deg\n', min(azimuth), max(azimuth));
    fprintf('Plage élévation: [%.2f, %.2f] deg\n', min(elevation), max(elevation));
    fprintf('Plage distance: [%.2f, %.2f] km\n', min(distance), max(distance));
    
    % Extraction des données pour la cible n
    idx = target_id == n;
    
    if ~any(idx)
        error('Cible %d non trouvée dans les données', n);
    end
    
    % Construction de la structure Target_data
    Target_data = struct();
    Target_data.id = n;
    Target_data.time = time(idx);
    Target_data.datetime = datetime_vals(idx);
    Target_data.x_true = x_true(idx);
    Target_data.y_true = y_true(idx);
    Target_data.z_true = z_true(idx);
    Target_data.azimuth = azimuth(idx);
    Target_data.elevation = elevation(idx);
    Target_data.distance = distance(idx);
    Target_data.n_measurements = sum(idx);
    
    % Traitement des mesures (si la fonction existe)
    if exist('Traitement_mesures', 'file')
        meas = Traitement_mesures(Target_data);
        
        Target_data.measurements = [meas.x_true; meas.y_true; meas.z_true];
        ECE=ecef2eci(utc,Target_data.measurements);
    end
    
    fprintf('\n=== CIBLE %d EXTRAITE ===\n', n);
    fprintf('Nombre de mesures: %d\n', Target_data.n_measurements);
end

