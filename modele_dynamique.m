function dot_etat = modele_dynamique(t, etat, sat_params)
    % modeled dynamique avec perturbations
    % etat: [x; y; z; vx; vy; vz] en km et km/s
    % sat_params: structure avec mu, J2, etc.
    
    r_vec = etat(1:3);
    v_vec = etat(4:6);
    r = norm(r_vec);
    
    %Gravité deux corps
    a_grav = -(sat_params.mu / r^3) * r_vec;
    
    % Calcul des perturbations
    a_pert = zeros(3,1);
    
    %J2
    if isfield(sat_params, 'J2') && sat_params.J2 ~= 0
        x = r_vec(1); y = r_vec(2); z = r_vec(3);
        R = sat_params.R_terre;
        coeff_J2 = (3/2) * sat_params.J2 * sat_params.mu * R^2 / r^5;
        a_J2 = coeff_J2 * [x*(5*z^2/r^2 - 1);
                           y*(5*z^2/r^2 - 1);
                           z*(5*z^2/r^2 - 3)];
        a_pert = a_pert + a_J2;
    end
    
    %Traînée atmosphérique
    if isfield(sat_params, 'BC') && sat_params.BC ~= 0
        % Calcul de l'altitude
        alt = r - sat_params.R_terre;
        
        % Modèle de densité exponentiel
        rho0 = 1.225e9;  % kg/km³ (équivalent 1.225 kg/m³)
        H = 8.5;         % km
        rho = rho0 * exp(-alt / H);
        
        % Vitesse relative (atmosphère co-rotative)
        omega_terre = [0; 0; 7.2921150e-5]; % rad/s
        v_atm = cross(omega_terre, r_vec);
        v_rel_vec = v_vec - v_atm;
        v_rel = norm(v_rel_vec);
        
        if v_rel > 0
            a_drag = -0.5 * sat_params.BC * rho * v_rel * v_rel_vec;
            a_pert = a_pert + a_drag;
        end
    end
    
    % Radiation solaire
    if isfield(sat_params, 'Cr') && sat_params.Cr ~= 0
        % Direction simplifiée du Soleil (pour test)
        u_soleil = [1; 0; 0];
        
        % Constante de radiation
        P_sol = 4.56e-6; % N/m²
        K_srp = P_sol * sat_params.Cr * (sat_params.A_srp / sat_params.masse) * 1e9;
        
        a_srp = -K_srp * u_soleil;
        a_pert = a_pert + a_srp;
    end
    
    % 3. Accélération totale
    a_tot = a_grav + a_pert;
    
    dot_etat = [v_vec; a_tot];
end
