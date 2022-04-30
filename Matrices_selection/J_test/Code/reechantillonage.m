function[sample]=reechantillonage(data)
% Heteroskedastic bootstrap (see hagemann, Journal of Econometrics 2012)
n=size(data,1);
e=randraw('rademacher', [], n);
de=diag(e);

sample=de*data;
% sample=zeros(n,1);
% choose=round(((n-1)*rand(n,1))+1);
% for j=1:n
% sample(j,:)=data(choose(j));
% end;


