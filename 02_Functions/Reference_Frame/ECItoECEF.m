function X_ecef = ECItoECEF(X_eci, TimeT, t0_datetime)
    % X_eci       : Matrice 6xN [x;y;z;vx;vy;vz] en ECI
    % TimeT       : Vecteur 1xN des secondes écoulées (0, 0.8, 1.6...)
    % t0_datetime : Objet datetime (ex: meas.datetime(1))
    
    N = size(X_eci, 2);
    X_ecef = zeros(6, N);
    
    for k = 1:N
        % On crée la date précise pour l'échantillon k
        current_dt = t0_datetime + seconds(TimeT(k));
        
        % On convertit en vecteur [AAAA, MM, JJ, HH, MM, SS] requis par la toolbox
        utc_vec = datevec(current_dt);
        
        r_eci_k = X_eci(1:3, k);
        v_eci_k = X_eci(4:6, k);
        
        % Transformation ECI -> ECEF
        % Cette fonction utilise en interne la rotation terrestre à cette date
        [r_ecef_k, v_ecef_k] = eci2ecef(utc_vec, r_eci_k, v_eci_k);
        
        X_ecef(1:3, k) = r_ecef_k;
        X_ecef(4:6, k) = v_ecef_k;
    end
end