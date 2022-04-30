function D = disteuc(x,y)

% Objectif : calcul la distance euclidienne entre chaque localisation
%% x : vecteur des longitudes
%% y : vecteur des lattitudes


n = length(x) ;   % nombre de localisations
T1=repmat(x',n,1);   % matrice où chaque colonne est identiquement composée des éléments de x
T2=repmat(y',n,1);   % matrice où chaque colonne est identiquement composée des éléments de y
T3=repmat(x,1,n);   % matrice où chaque ligne est identiquement composée des éléments de x
T4=repmat(y,1,n);   % matrice où chaque ligne est identiquement composée des éléments de y
D=sqrt((T3-T1).^2+(T4-T2).^2);   % calcul de la distance euclidienne