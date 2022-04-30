function W = faire_voisin(xc,yc,m)

% Meme programme que Lesage sauf qu'il utilise trouve_voisin (avec distance
% arc).
%% xc = vecteur des lattitudes
%% yc = vecteur des longitudes
%% m = nombre de voisins


if nargin == 3
[n junk] = size(xc);    
else,
error('make_neighborsw: Wrong # of input arguments');
end;


nnlist = trouve_voisin(xc,yc,m);

% convert the list into a row-standardized spatial weight matrix
rowseqs=(1:n)';
vals1=ones(n,1)*(1/m);
vals0=zeros(n,1);

for i=1:m;

colseqs=nnlist(:,i);
ind_to_keep=logical(colseqs>0);

z1=[rowseqs colseqs vals1];
z1=z1(ind_to_keep,:);

z2=[rowseqs rowseqs vals0];
%this last statement makes sure the dimensions are right
z=[z1
   z2];

if i == 1
    W = spconvert(z);
else
    W = W + spconvert(z);
end;

end;


