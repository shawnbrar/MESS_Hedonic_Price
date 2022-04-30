function Wd = dinvarc(x,y,k)

% Objectif : calcul de la matrice de distance arc inverse standardisée
% entre chaque localisation
%% x : vecteur des latitudes
%% y : vecteur des longitudes
%% k : puissance

n = length(x);
D = distarc(x,y);% calcul de la distance arc entre chaque localisation
% D=D/1000; % dividing each element by 1000
E = D + eye(n);   % met des 1 sur la diagonale pour éviter de diviser par 0
F = (ones(n)./E).^k;   % calcul les distances inverses
Wd= F - eye(n);   % enlève les 1 sur la diagonale
% Wd = normw(G);   % normalise la matrice