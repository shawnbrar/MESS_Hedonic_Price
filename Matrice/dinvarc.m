function Wd = dinvarc(x,y,k)

% Objectif : calcul de la matrice de distance arc inverse standardisée
% entre chaque localisation
%% x : vecteur des longitudes
%% y : vecteur des lattitudes
%% k : puissance

n = length(x);
D = distarc(x,y);   % calcul de la distance arc entre chaque localisation
E = D + eye(n);   % met des 1 sur la diagonale pour éviter de diviser par 0
F = (ones(n)./E).^k;   % calcul les distances inverses
Wd= F - eye(n);   % enlève les 1 sur la diagonale
% Wd = normw(G);   % normalise la matrice