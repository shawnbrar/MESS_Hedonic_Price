function Wd = dexparc(x,y,k)

% Objectif : calcul de la matrice de distance arc exponentielle inverse standardisée
% entre chaque localisation
%% x : vecteur des longitudes
%% y : vecteur des lattitudes
%% k : puissance

n = length(x);
D = distarc(x,y);    % calcul de la distance arc entre chaque localisation
E = D / 1000;   % met des 1 sur la diagonale pour éviter de diviser par 0
F = exp(-k*E);
Wd = F - eye(n);   % enlève les 1 sur la diagonale
% Wd = normw(G);   % normalise la matrice