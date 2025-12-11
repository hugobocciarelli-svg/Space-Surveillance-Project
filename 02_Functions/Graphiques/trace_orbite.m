function trace_orbite(Etats)

    x = Etats(1, :);
    y = Etats(2, :);
    z= Etats(3, :);
    R_T = 6378; % km

    % Création de la sphère Terre
    [Xe, Ye, Ze] = sphere(50);
    figure;
    hold on;
    grid on;

    %Terre
    surf(Xe * R_T, Ye * R_T, Ze * R_T, 'FaceColor',[0.2 0.4 1], ...
         'EdgeColor','none','FaceAlpha',0.7);

    % Orbite du satellite
    plot3(x, y, z, 'r', 'LineWidth', 2);

    % Vue 
    axis equal;
    xlabel('x (km)');
    ylabel('y (km)');
    zlabel('z (km)');
    title('Trajectoire orbitale autour de la Terre');
    
    view(45, 25);
    hold off;

end