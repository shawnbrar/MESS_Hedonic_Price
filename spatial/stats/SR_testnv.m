function SRtest = SR_testnv(y,W,nv1,nv2,nv3,q,cont)
%%Computes the spatial runs statistic
% May 2015
% This is the last version of the code. We need to disregard the SR_test.m
% y is the continuous data set. If the data is categorical then coment lines 21 to 35
% and line 44.
% W is the association squeme.
% q is the number of categories.
% nvar number of types of covariances for the association esqueme W
% lat latitude coordinates for W
% long longitude coordinates for W

n=length(y);

% Cont is a binary variable that takes on the value of 1 if data are
% continuous and 0 if data are categorical.

m=zeros(q,1);
pprod=zeros(q,1);

if cont==1
for qq=1:q-1;
    cuts(qq)=quantile(y,qq/q);
end;
%mex=median(xx);
yy{1}=y<=cuts(1);
if q >2
for j=2:q-1
yy{j}=j*(cuts(j-1)<y & y<=cuts(j));
end
end
yy{q}=q*(y>cuts(q-1));
% for j=1:q-2;
%     yy{j+1}=(j+1)*(cuts(j)<=xx & xx<=cuts(j+1));
% end;

y=zeros(n,1);

for k=1:q;
%         y(find(yy{k}))=k+zeros(length(find(yy{k})),1)
        y(find(yy{k}))=k;
        m(k)=length(find(yy{k}));
        pprod(k)=m(k)*(n-m(k));
    end;
else 
    
    for k=1:q;
%         y(find(yy{k}))=k+zeros(length(find(yy{k})),1)
%         y(find(yy{k}))=k;
        m(k)=length(find(y==k));
        pprod(k)=m(k)*(n-m(k));
    end;
end


% here we categorize the original data set y into the q categories
%compute the m_k needed for the computation of mean and variance
%pprod is needed for the computation of p

p=sum(pprod)/(n*(n-1));

    
%here we find out the neighbors and sort them for each location
for i=1:n;
    nn{i}=find(W(i,:))';
    ns(i)=length(nn{i});
    ww{i}=W(i,nn{i});
    %dis{i}=(lat(i)-lat(nn{i})).^2+(long(i)-long(nn{i})).^2; % Sort according to
%    Euclidean distance
    [di aa]=sort(-ww{i}); % Sort according to the biggest weights (as weights are all positive, we take minus so taht the biggest
%     positive will become thee smallest negative and the sort function ca
%     be used without modification. 
%   [di aa]=sort(dis{i});
    nb{i}=nn{i}(aa);
end;

%Here we compute the runs starting at each location and it sum is the total number of runs 
for i=1:n;
%     ll{i}=find(y(nb{i})==y(i));
%     ru=ones(length(nb{i}),1);
%     ru(ll{i})=0;
%     run{i}=ru;
    rrun{i}=abs(diff(y(vertcat(i,nb{i}))));
    runs(i)=1+length(find(rrun{i}>0));
    %syi{i}=
end;
SR=sum(runs);
%The mean of the statistic
meanSR=n+p*sum(ns);


%%%%%COMPUTING THE VARIANCE%%%%%
%% case 1 %%
aux1=zeros(q,1);
aux2=zeros(q,1);
aux31=zeros(q,1);
for k=1:q
%     aux1(k)=m(k)*(n-m(k));
    aux2(k)=m(k)*(n-m(k))*(n-m(k)-1);
    aux31(k)=m(k)*(m(k)-1)*(n-m(k))*(n-m(k)-1);
end

% sum1=sum(aux1);
sum1=sum(pprod);
sum2=sum(aux2);
sum3=2*sum(aux31);
% t2=0;
% aux3=zeros(q^2,1);

  t3=0;
  %aux13=zeros(q^3,1);
  aux32=zeros(q^3,1);
    for k=1:q;
     for c=1:q;
       for d=1:q;
           t3=t3+1;
          aux32(t3)=m(k)*m(c)*m(d)*(n-m(d)-2);
          %aux33(t3)=m(k)*(m(k)-1)*m(c)*m(d)*(m(d)-1);
          
         if c==k
             aux32(t3)=0;
%              aux33(t3)=0;
             
         end; 
         if d==k
             aux32(t3)=0;
%              aux33(t3)=0;
             
         end;
         if d==c
             aux32(t3)=0;
%              aux33(t3)=0;
         end;
       end;
          end;
 
   end;
 sum32=sum(aux32);
%  sum33=sum(aux33);
%  sum42=2*sum33;
%  
 
 
            
  var1=1/(n*(n-1)*(n-2)*(n-3))*(sum3+sum32);
  var2=1/(n*(n-1)*(n-2))*(sum2);
  var3=1/(n*(n-1))*(sum1);
%   var4=1/(n*(n-1)*(n-2)*(n-3)*(n-4))*(sum41+sum42);
%   var5=1/(n*(n-1)*(n-2)*(n-3))*(sum51);
%   var6=1/(n*(n-1)*(n-2)*(n-3))*(sum61);
%   var7=2/(n*(n-1)*(n-2)*(n-3))*(sum71);
%   var8=1/(n*(n-1)*(n-2))*(sum81);
  

    
  varSR=p*(1-p)*sum(ns)+2*sum(ns-1)*(var2-p^2)+sum((ns-1).*(ns-2))*(var1-p^2)+...
  nv1*var1+nv2*var2+nv3*var3-(nv1+nv2+nv3)*p^2;
  
  %the test statistic which is N(0,1) distributed
  SRtest=(SR-meanSR)/sqrt(varSR);
  
%   result.m=m';
%   result.mean=meanSR;
%   result.var=varSR;
