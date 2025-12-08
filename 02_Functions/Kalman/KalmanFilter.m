function [x_updated,P_updated] = KalmanFilter(x_pred,P_pred,Z,H,V)
% Kaklman filter : La fonction prend en argument les prédictions du vecteurs
% d'état [X]= [x,y,z,vx,vy,vz] dans ECI, celle de la matrice de covariance [P] , la mesure ainsi que la
% matrice [H] (Matrice des mesures) ainsi que [V] le bruit associé à la mesure 
% Les sorties de la fonction est le l'etat X mise à jour ainsi que la
% matrice P mise à jour aussi.
% Les équations d'etats sont décrite tel que : {Xk+1=Fk * Xk + Wk et Zk= HkZk
% +Vk)

% Correction du ain de Kalman
K = P_pred * H'* inv(H * P_pred * H' + V);

% Mise à jour
x_updated = x_pred + K * (Z-H*x_pred);
P_updated = (eye(6) - K * H) * P_pred;
end
