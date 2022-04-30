function D = distarc(x,y)

% Objectif : calcul la distance sph�rique entre chaque localisation
%% x : vecteur des latitudes
%% y : vecteur des longitudes
%  For further info concerning lat-long conversion, see the following site
%  http://www.uwgb.edu/dutchs/usefuldata/utmformulas.htm
% lat and long are assumed to be in decimal degree and not in degrees
% minutes seconds
% See http://geodesie.ign.fr/contenu/fichiers/Distance_longitude_latitude.pdf for further details

n = length(x) ;   % nombre de localisations
% Transformation en radian
X = (pi/180)*x; 
Y = (pi/180)*y;
T1=repmat(X',n,1);   % matrice o� chaque colonne est identiquement compos�e des �l�ments de x
T2=repmat(Y',n,1);   % matrice o� chaque colonne est identiquement compos�e des �l�ments de y
T3=repmat(X,1,n);   % matrice o� chaque ligne est identiquement compos�e des �l�ments de x
T4=repmat(Y,1,n);   % matrice o� chaque ligne est identiquement compos�e des �l�ments de y
% Distance angulaire en gradiant S_a,b = arc cos(sin(lat_a)*sin(lat_b) + cos(lat_a)*cos(lat_b)*cos(long_b-long_a))
% Pour convertir cette distance en mettre, on multiplie S_a,b par un rayon
% de la terre conventionnel (exprim� en m�tres): 6 378 137 m) 
E = 6378 * acos(cos(abs(T2-T4)).*cos(T1).*cos(T3) + sin(T1).*sin(T3));
% Elimination des �l�ments de la diagonale
D = E.*(ones(n,n)-eye(n,n));
% La distance est calcul�e en km car le rayon est exprim� en km.