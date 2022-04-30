%Project: Contagion in OECD's government bond market
%Purpose: This function apply a specific normalization to a weight
%         matrix. We consider three kind of normalization:
%         - row-normalization (meth=0)
%         - Normalization by the max eigenvalue (meth=1)
%         - Normalization by the max sum of the rows and columns (meth=2)
%
%Input: - W: origianl weigth matrix
%       - meth: scalar 
%
%Output: Normalized weigth matrix
%           
%Author:  Cyrille D.
%created: 26/06/2014
%Last update: 06/08a/2014

function [Ws]=Normaliz(W,meth)
    
    if meth==0                      %Row-Normalization
        rsum=repmat(sum(W,2),1,size(W,2));
        rsum(rsum==0)=1;
        Ws=W./rsum;
    elseif meth==1                  %Normalization by the max eigenvalue
        tau=max(eig(W));
        Ws=(1/tau)*W;
    elseif meth==2                  %Normalization by the max sum of the rows and columns
        sumi=sum(abs(W),2);
        sumj=sum(abs(W),1);
        tau_star=min(max(sumi),max(sumj));
        Ws=(1/tau_star)*W;
    end
    
end
 
