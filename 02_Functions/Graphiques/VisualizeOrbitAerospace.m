function VisualizeOrbitAerospace(State, Time_list,meas,t_sim)   
    rv_eci_6xN = State;      
    t_sec      = Time_list;   
    rv_eci = rv_eci_6xN';           % N × 6
    t_relatif = Time_list - Time_list(1);
    
    r_eci = rv_eci(:, 1:3);             % N × 3 km
    v_eci = rv_eci(:, 4:6);             % N × 3 km/s
    
    % Conversion unités
    r_eci_m = r_eci * 1000;             % N × 3 mètres
    v_eci_ms = v_eci * 1000;            % N × 3 m/s
    
    % Création timeseries 
    position_ts = timeseries(r_eci_m, t_relatif, 'Name', 'Position ECI');
    velocity_ts = timeseries(v_eci_ms, t_relatif, 'Name', 'Velocity ECI');

    % Interpolation linéaire 
    position_ts.DataInfo.Interpolation = "linear";  
    velocity_ts.DataInfo.Interpolation = "linear";

    % Scénario 
    startTime = meas.datetime(1);  
    sc = satelliteScenario(startTime, startTime + seconds(t_sim), 30);
    
    % Satellite
    sat = satellite(sc, position_ts, velocity_ts,"Name", "Trajectoire ECI Custom");
    sat.MarkerColor = [1 0.8 0];   
    sat.MarkerSize  = 10;
    
    % Lancement du viewer 
    play(sc);
end