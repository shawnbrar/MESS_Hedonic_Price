function f=wrook(n)
nn=n*n;
W=zeros(nn,nn);
for i=1:nn
    for j=1:nn
        if abs(i-j)==1|abs(i-j)==n
            W(i,j)=1;
        end
    end
end

tempw=0;
for i=1:nn;
    for j=1:nn;
        tempw=tempw+W(i,j);
    end
    W(i,:)=W(i,:)/tempw; %W is normalized 
    tempw=0;
end

f=W;



            