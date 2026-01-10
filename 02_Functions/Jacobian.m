
function J = Jacobian(state)
mu = 398600.4418;  
J2 = 0.0010836; 
Re = 6371;



%Calcule la matrice Jacobienne 6x6 du modèle orbital avec J2
x  = state(1);
y  = state(2);
z  = state(3);

%Calcul préalable des puissances de r pour limiter le temps de calcul
r2 = x^2 + y^2 + z^2;
r  = sqrt(r2);
r3 = r^3;
r5 = r^5;
r7 = r^7;

% On pose une constante C
C = 1.5 * mu * J2 * Re^2;   

% Initialisation 
J = zeros(6,6);
J(1,4) = 1;
J(2,5) = 1;
J(3,6) = 1;

% ∂ax/∂x 
J(4,1) = (-mu/r3 + 3*mu*x^2/r5) ...
    + (C/r5) * (-1 + 5*(z^2/r2) + 5*(x^2/r2)*(1 - 7*(z^2/r2)));

% ∂ay/∂y  
J(5,2) = (-mu/r3 + 3*mu*y^2/r5) ...
    + (C/r5) * (-1 + 5*(z^2/r2) + 5*(y^2/r2)*(1 - 7*(z^2/r2)));

% ∂az/∂z 
J(6,3) = (-mu/r3 + 3*mu*z^2/r5) ...
    + (C/r5) * (-3 + 30*(z^2/r2) - 35*(z^2/r2)^2);

% ∂ax/∂y = ∂ay/∂x
J(4,2) = 3*mu*x*y/r5 + 5*C*x*y/r7 * (1 - 7*(z^2/r2));
J(5,1) = J(4,2);

% ∂ax/∂z = ∂az/∂x
J(4,3) = 3*mu*x*z/r5 + 5*C*x*z/r7 * (3 - 7*(z^2/r2));
J(6,1) = J(4,3);

% ∂ay/∂z = ∂az/∂y
J(5,3) = 3*mu*y*z/r5 + 5*C*y*z/r7 * (3 - 7*(z^2/r2));
J(6,2) = J(5,3);


end