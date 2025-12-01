%% SETUP      
t0 = 0;
dt = 10;
t_end = 6000;
TimeT = t0:dt:t_end;
n_steps = length(TimeT);

%% Paramètres initiaux
X0=[7078;0;0;0;7.5;0];

% Initialisation du vecteur d'état
Etats = zeros(6, n_steps);
Etats(:, 1)=X0;


for i = 1:n_steps-1
    % Mise à jour
    t = TimeT(i); % Update time for the current step
    Etats(:, i+1) = RungeKutta4(@modele_dynamique,t,Etats(:, i),dt);
end

%Tracer l'orbite
trace_orbite(Etats);