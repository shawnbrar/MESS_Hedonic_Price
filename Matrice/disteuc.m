function D = disteuc(x,y)

% Objectif : calcul la distance euclidienne entre chaque localisation
%% x : vecteur des longitudes
%% y : vecteur des lattitudes


n = length(x) ;   % nombre de localisations
T1=repmat(x',n,1);   % matrice o� chaque colonne est identiquement compos�e des �l�ments de x
T2=repmat(y',n,1);   % matrice o� chaque colonne est identiquement compos�e des �l�ments de y
T3=repmat(x,1,n);   % matrice o� chaque ligne est identiquement compos�e des �l�ments de x
T4=repmat(y,1,n);   % matrice o� chaque ligne est identiquement compos�e des �l�ments de y
D=sqrt((T3-T1).^2+(T4-T2).^2);   % calcul de la distance euclidienne