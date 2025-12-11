
function PlotTrajectoryComparison(TimeT, Corrected_states, Position, Measures,meas_time, meas_datetime)
    time_meas_col = meas_datetime(:);
    
    % 1. Calculer les positions dans ECI 
    Position_ECI = ECEFtoECI(Position, time_meas_col);
    
    % 2. Extraire la position en sortie du filtre
    Calc_ECI = Corrected_states(1:3, :);
    

    figure;
    sgtitle('Comparaison des Trajectoires Calculée et Mesurée (Repère ECI)');
    Composantes = {'X (km)', 'Y (km)', 'Z (km)'};
    
    for k = 1:3
        subplot(3, 1, k);
        plot(TimeT, Calc_ECI(k, :), 'b-', 'LineWidth', 1.5);

        hold on;
        plot(meas_time, Position_ECI(k, :), 'ro', 'MarkerSize', 2, 'MarkerFaceColor', 'r');
        plot(meas_time, Measures(k, :), 'ro', 'MarkerSize', 2, 'MarkerFaceColor', 'g');
        grid on;
        ylabel(Composantes{k});
        if k == 1
            legend('Trajectoire Calculée', 'Mesures réelles', 'Location', 'best');
        end
        if k == 3
            xlabel('Temps (s)');
        end
        
        hold off;
    end
end