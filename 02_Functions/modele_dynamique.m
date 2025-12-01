function dot_etat = modele_dynamique(t,etat)
%Two body problem model
    mu =398600.4418; %km^3/s^2
    r_vec = etat(1:3); % Vecteur position [x; y; z]
    v_vec = etat(4:6); % Vecteur vitesse [vx; vy; vz]
    r = norm(r_vec);
    a_gravite = -(mu / r^3) * r_vec;
   
    dot_etat = [
        v_vec;    
        a_gravite 
    ];
end
