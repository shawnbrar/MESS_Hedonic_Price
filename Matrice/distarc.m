function D = distarc(x,y)

% Objectif : calcul la distance sph�rique entre chaque localisation
%% x : vecteur des longitudes
%% y : vecteur des lattitudes

n = length(x) ;   % nombre de localisations
% Transformation en radian
X = (pi/180)*x;
Y = (pi/180)*y;
T1=repmat(X',n,1);   % matrice o� chaque colonne est identiquement compos�e des �l�ments de x
T2=repmat(Y',n,1);   % matrice o� chaque colonne est identiquement compos�e des �l�ments de y
T3=repmat(X,1,n);   % matrice o� chaque ligne est identiquement compos�e des �l�ments de x
T4=repmat(Y,1,n);   % matrice o� chaque ligne est identiquement compos�e des �l�ments de y
E = 6378 * acos(cos(abs(T2-T4)).*cos(T1).*cos(T3) + sin(T1).*sin(T3));
% Elimination des �l�ments de la diagonale
D = E.*(ones(n,n)-eye(n,n));