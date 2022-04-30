function [nvar1, nvar2, nvar3]=creation_nvar_SR(W)

nvarf1=0;
nvarf2=0;
nvarf3=0;
n=size(W,1);
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
     nnb{i}=vertcat(i,nb{i});
end;

pq=0;
tic;
%parpool(7);
for i=1:n
   
    
%     for i = 1:n
    for j=1:n
        if j==i
                 
        else
                 pq=pq+1 ;  
          
        for k=2:ns(i)+1
          csi=[nnb{i}(k-1), nnb{i}(k)];
          for c=2:ns(j)+1
           ctj=[nnb{j}(c-1), nnb{j}(c)];
              
           if  isempty(intersect(csi,ctj))==1
                  
              nvarf1=nvarf1+1;
           end;
              
              if length(intersect(csi,ctj))==1
                 nvarf2=nvarf2+1;   
              end;
              if length(intersect(csi,ctj))==2
                  nvarf3=nvarf3+1;
              end
              
              
              end;
              
             
              
          end;
      end;
    end;
end;
nvar1=nvarf1;
nvar2=nvarf2;
nvar3=nvarf3;
toc;

%delete(gcp);
