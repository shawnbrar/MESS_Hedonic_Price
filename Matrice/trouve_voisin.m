function nnlist = trouve_voisin(xc,yc,m)

%% xc = vecteur des longitudes
%% yc = vecteur des lattitudes



if nargin ~= 3
error('find_neighbors: 3 input arguments required');  % Vérification du nombre d'arguments
end;

% Vérification des erreurs
[n junk] = size(xc);
if junk ~= 1
xc = xc';
end;
[n2 junk] = size(yc);
if junk ~= 1
yc = yc';
end;
if n ~= n2
error('find_neighbors: xc,yc inputs must be same size');
end;

nnlist = zeros(n,m);   % Création d'un matrice de 0 qu'il faut remplir

for i=1:n;
n = length(xc) ;   % nombre de localisations
% Transformation en radian
X = (pi/180)*xc;
Y = (pi/180)*yc;
    xi = X(i,1);
    yi = Y(i,1);
E = 6378 * acos(cos(abs(yi-Y)).*cos(xi).*cos(X) + sin(xi).*sin(X)); % Calcul la distance Arc

[xds xind] = sort(E);        % tri de manière croissante avec la distance et renvoie deux matrices :
                                  % une avec les valeurs (xds), l'autre avec les
                                  % numéros de lignes (xind).
nnlist(i,1:m) = xind(2:m+1,1)';   % Pour la ième colonne de nnlist, on prend les 10 premiers (sauf le 0)
end;

