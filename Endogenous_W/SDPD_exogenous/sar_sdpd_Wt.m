function results = sar_sdpd_Wt(y,x,W,info)
timet = clock; % start the clock for overall timing
% default options
options = optimset('fminsearch');
parm = 0;

fields = fieldnames(info);
nf = length(fields);
for i=1:nf
    if strcmp(fields{i},'parm')
        parm = info.parm;
    elseif strcmp(fields{i},'ted')
        ted = info.ted;
    elseif strcmp(fields{i},'tl')
        tl = info.tl; % tar lag index
    elseif strcmp(fields{i},'stl')
        stl = info.stl; % star lag index
    end;
end;

n=info.n;
[nt,kx]=size(x);
t=nt/n;
y2=y;
x2=x;

if ted==1
    Jn=speye(n)-1/n*ones(n,1)*ones(1,n);
    [Fnn junk]=eig(Jn);
    F=Fnn(:,2:n);
    
    y2=reshape(y2,n,t+1);y2temp=zeros(n-1,t+1);
    y2temp=F'*y2;
    y2=reshape(y2temp,(n-1)*(t+1),1);
    
    if isempty(x)==0
        x2=reshape(x2,n,t,kx);x2temp=zeros(n-1,t,kx);
        for i=1:kx
            x2temp(:,:,i)=F'*x2(:,:,i);
        end
        x2=reshape(x2temp,(n-1)*t,kx);
    else
        x2=[];
    end
    
    W2=zeros(n-1,(n-1)*(t+1));
    for i=1:t+1
        W2(:,1+(i-1)*(n-1):i*(n-1))=F'*W(:,1+(i-1)*n:i*n)*F;
    end
    W=W2;
    n=n-1;
    nt=n*t;
end

%%newly added for ted=2
if ted==2
    y2temp=y2;
    x2temp=x2;

    Jn=speye(n)-1/n*ones(n,1)*ones(1,n);
    for i=1:t+1
        y2temp(1+(i-1)*n:i*n,:)=Jn*y2(1+(i-1)*n:i*n,:);

    end
    for i=1:t
        x2temp(1+(i-1)*n:i*n,:)=Jn*x2(1+(i-1)*n:i*n,:);

    end
    y2=y2temp;
    x2=x2temp;

end

yt=y2(n+1:n+nt);
ytl=y2(1:nt);
ysl=zeros(nt,1);ystl=zeros(nt,1);
for i=1:t
    ysl(1+(i-1)*n:i*n)=W(:,1+i*n:(i+1)*n)*yt(1+(i-1)*n:i*n);
    ystl(1+(i-1)*n:i*n)=W(:,1+(i-1)*n:i*n)*ytl(1+(i-1)*n:i*n);
end

yt=reshape(yt,n,t);
temp=mean(yt')';
temp=temp*ones(1,t);
yt=yt-temp;
yt=reshape(yt,nt,1);

ytl=reshape(ytl,n,t);
temp=mean(ytl')';
temp=temp*ones(1,t);
ytl=ytl-temp;
ytl=reshape(ytl,nt,1);

ysl=reshape(ysl,n,t);
temp=mean(ysl')';
temp=temp*ones(1,t);
ysl=ysl-temp;
ysl=reshape(ysl,nt,1);

ystl=reshape(ystl,n,t);
temp=mean(ystl')';
temp=temp*ones(1,t);
ystl=ystl-temp;
ystl=reshape(ystl,nt,1);

xt=x2;
if isempty(x) == 0
    [junk,kx]=size(xt);
    xt=reshape(xt,n,t,kx);
    for i=1:kx
        temp=mean(xt(:,:,i)')';
        temp=temp*ones(1,t);
        xt(:,:,i)=xt(:,:,i)-temp;
    end
    xt=reshape(xt,nt,kx);
else
    xt=[];
end

if stl + tl == 2
    zt=[ytl ystl xt];
elseif stl + tl == 1
    if stl == 1, zt=[ystl xt]; else zt=[ytl xt]; end
elseif stl + tl == 0
    error('Wrong Info input,Our model has dynamic term anyway');
else
    error('Double-Check stl & tl # in Info structure ');
end

[junk kz]=size(zt);
f_sar_sdpd_Wt_oc=@(x) f_sar_sdpd_Wt(x,yt,ysl,zt,W);
[pout,like,exitflag,output]=fminsearch(f_sar_sdpd_Wt_oc,parm,options);

% if exitflag == 0
% fprintf(1,'\n sac: convergence not obtained in %4d iterations \n',output.iterations);
% end;
results.iter = output.iterations;

results.lik = -t*like;% see f_sar_sdpd_Wt

rho = pout(1,1);
rho = 2/(1+exp(-rho))-1;

Syt=yt-rho*ysl;
b0 = zt\Syt;
e = Syt- zt*b0;
results.beta = b0;
results.rho = rho;
results.resid = e;
results.yhat = yt-e;
sigu = e'*e;
sige = sigu/nt;
results.sige = sige;
results.theta=[results.beta;results.rho;results.sige];
theta=results.theta;

[junk nvar]=size(zt);
bhat = results.beta;

xpx = zeros(nvar+2,nvar+2);
xpx(1:nvar,1:nvar) = (1/sige)*(zt'*zt);
trG=0;
trG2=0;
for i=1:t
    St=speye(n)-rho*W(:,1+i*n:(i+1)*n);
    Wt=W(:,1+i*n:(i+1)*n);
    Gt=Wt*inv(St);
    trG=trG+trace(Gt);
    trG2=trG2+trace(Gt^2);
end

term1 = trG2;
term2 = (1/sige)*ysl'*ysl;
xpx(nvar+1,nvar+1) = term1+term2;
xpx(nvar+2,nvar+2) = nt/(2*sige*sige);
xpx(1:nvar,nvar+1) = (1/sige)*(zt'*ysl);
xpx(nvar+1,1:nvar) = xpx(1:nvar,nvar+1)';
xpx(nvar+1,nvar+2) = (1/sige)*trG;
xpx(nvar+2,nvar+1) = xpx(nvar+1,nvar+2)';

xpx=xpx/nt;
xpxi = invpd(xpx);
tmp = diag(abs(xpxi));

bvec = results.theta;
results.tstat = sqrt(n*t)*bvec./sqrt(tmp(1:nvar+2,1));

results.std=sqrt(tmp)/sqrt(n*t);
results.y = y;

In=eye(n);

bias=zeros(kz+2,1);
S=kron(ones(1,t+1),In)-rho*W;%first S is time period 0
G=zeros(n,n*t);
for s=1:t
    G(:,1+(s-1)*n:s*n)=W(:,1+s*n:(s+1)*n)*inv(S(:,1+s*n:(s+1)*n));
end

A=zeros(n,n*t);%first A is time period 1

sumG=zeros(n,n);
for s=1:t
    sumG=sumG+ G(:,1+(s-1)*n:s*n);
end

for s=1:t
    if stl + tl == 2
        gamma_coff=b0(1);rho_coff=b0(2);
        A(:,1+(s-1)*n:s*n)=inv(S(:,1+s*n:(s+1)*n))*(gamma_coff*In+rho_coff*W(:,1+(s-1)*n:s*n));
    elseif stl + tl == 1
        if stl == 1
            rho_coff==b0(1);
            A(:,1+(s-1)*n:s*n)=inv(S(:,1+s*n:(s+1)*n))*(rho_coff*W(:,1+(s-1)*n:s*n));
        else
            gamma_coff=b0(1);
            A(:,1+(s-1)*n:s*n)=inv(S(:,1+s*n:(s+1)*n))*(gamma_coff*In);
        end
    elseif stl + tl == 0
        error('Wrong Info input,Our model has dynamic term anyway');
    else
        error('Double-Check stl & tl # in Info structure ');
    end
end

bias1=0;bias2=0;bias3=0;
if stl + tl == 2
    for s=1:t-1
        temp1a=eye(n);%this is actually h=1
        temp2a=W(:,1+s*n:(s+1)*n);
        temp3a=G(:,1+s*n:(s+1)*n)*(gamma_coff*temp1a+rho_coff*temp2a);
        for h=2:t-s
            temp1b=eye(n);
            for j=1:h-1
                temp1b=temp1b*A(:,1+(s+j-1)*n:(s+j)*n);
            end
            temp1a=temp1a+temp1b;
            temp2a=temp2a+W(:,1+(s+h-1)*n:(s+h)*n)*temp1b;
            temp3a=temp3a+gamma_coff*G(:,1+(s+h-1)*n:(s+h)*n)*temp1b+rho_coff*G(:,1+(s+h-1)*n:(s+h)*n)*W(:,1+(s+h-1)*n:(s+h)*n)*temp1b;
        end
        bias1=bias1+trace(inv(S(:,1+s*n:(s+1)*n))*temp1a);
        bias2=bias2+trace(inv(S(:,1+s*n:(s+1)*n))*temp2a);
        bias3=bias3+trace(inv(S(:,1+s*n:(s+1)*n))*temp3a);
    end
    
    bias(1,1)=(1/nt)*bias1;
    bias(2,1)=(1/nt)*bias2;
    bias(kz+1,1)=(1/nt)*bias3+(1/nt)*trG;
    bias(kz+2,1)=inv(2*sige);
    
    
elseif stl + tl == 1
    
    if stl == 1
        for s=1:t-1
            
            temp1a=eye(n);%this is actually h=1
            temp2a=W(:,1+s*n:(s+1)*n);
            temp3a=G(:,1+s*n:(s+1)*n)*rho_coff*temp2a;
            for h=2:t-s
                temp1b=eye(n);
                for j=1:h-1
                    temp1b=temp1b*A(:,1+(s+j-1)*n:(s+j)*n);
                end
                temp2a=temp2a+W(:,1+(s+h-1)*n:(s+h)*n)*temp1b;
                temp3a=temp3a+rho_coff*G(:,1+(s+h-1)*n:(s+h)*n)*W(:,1+(s+h-1)*n:(s+h)*n)*temp1b;
            end
            bias2=bias2+trace(inv(S(:,1+s*n:(s+1)*n))*temp2a);
            bias3=bias3+trace(inv(S(:,1+s*n:(s+1)*n))*temp3a);
        end
        
        bias(1,1)=(1/nt)*bias2;
        bias(kz+1,1)=(1/nt)*bias3+(1/nt)*trG;
        bias(kz+2,1)=inv(2*sige);
        
        
    else
        for s=1:t-1
            
            temp1a=eye(n);%this is actually h=1
            temp2a=W(:,1+s*n:(s+1)*n);
            temp3a=G(:,1+s*n:(s+1)*n)*gamma_coff*temp1a;
            for h=2:t-s
                temp1b=eye(n);
                for j=1:h-1
                    temp1b=temp1b*A(:,1+(s+j-1)*n:(s+j)*n);
                end
                temp1a=temp1a+temp1b;
                temp3a=temp3a+gamma_coff*G(:,1+(s+h-1)*n:(s+h)*n)*temp1b;
            end
            bias1=bias1+trace(inv(S(:,1+s*n:(s+1)*n))*temp1a);
            bias3=bias3+trace(inv(S(:,1+s*n:(s+1)*n))*temp3a);
        end
        
        bias(1,1)=(1/nt)*bias1;
        bias(kz+1,1)=(1/nt)*bias3+(1/nt)*trG;
        bias(kz+2,1)=inv(2*sige);
        
    end
elseif stl + tl == 0
    error('Wrong Info input,Our model has dynamic term anyway');
else
    error('Double-Check stl & tl # in Info structure ');
end

Bias2=zeros(kz+2,1);
Bias2(kz+1,1)=(1/nt)*ones(1,n)*sumG*ones(n,1);;
Bias2(kz+2,1)=inv(2*sige);

if ted==1
    theta1=theta+xpxi*bias/t;
end

if ted==2
    theta1=theta+xpxi*bias/t+xpxi*Bias2/n;
end

rho=theta1(nvar+1);
b0=theta1(1:nvar);
sige=theta1(nvar+2);

xpx = zeros(nvar+2,nvar+2);
xpx(1:nvar,1:nvar) = (1/sige)*(zt'*zt);
trG=0;
trG2=0;
for i=1:t
    St=speye(n)-rho*W(:,1+i*n:(i+1)*n);
    Wt=W(:,1+i*n:(i+1)*n);
    Gt=Wt*inv(St);
    trG=trG+trace(Gt);
    trG2=trG2+trace(Gt^2);
end

term1 = trG2;
term2 = (1/sige)*ysl'*ysl;
xpx(nvar+1,nvar+1) = term1+term2;
xpx(nvar+2,nvar+2) = nt/(2*sige*sige);
xpx(1:nvar,nvar+1) = (1/sige)*(zt'*ysl);
xpx(nvar+1,1:nvar) = xpx(1:nvar,nvar+1)';
xpx(nvar+1,nvar+2) = (1/sige)*trG;
xpx(nvar+2,nvar+1) = xpx(nvar+1,nvar+2)';

xpx=xpx/nt;
xpxi = invpd(xpx);
tmp = diag(abs(xpxi));

SIGi1=xpxi;
results.SIGi1=SIGi1;
SIG1=invpd(SIGi1);
results.SIG1=SIG1;

bvec = theta1;
results.tstat1 = sqrt(n*t)*bvec./sqrt(tmp(1:nvar+2,1));

std1=sqrt(tmp)/sqrt(n*t);

results.theta1=theta1;
results.std1=std1;


%%following is to compute the three impact--direct, total, and indirect

bhat=theta1(tl+stl+1:tl+stl+kx);

f_direct=zeros(kx,t);
f_total=zeros(kx,t);
f_indirect=zeros(kx,t);

sd_f_direct=zeros(kx,t);
sd_f_total=zeros(kx,t);
sd_f_indirect=zeros(kx,t);

for tempt=1:t
    for tempk=1:kx
          Wt=W(:,1+tempt*n:(tempt+1)*n);
          St=speye(n)-rho*Wt;
          Sti=inv(St);
          Rt= Sti*bhat(tempk);
          f_direct(tempk,tempt)=1/n*trace(Rt);
          f_total(tempk,tempt)=1/n*sum(sum(Rt));
          f_indirect(tempk,tempt)=1/n*trace(Rt);-1/n*sum(sum(Rt));
          
          P_f_direct=zeros(kz+2,1);
          P_f_total=zeros(kz+2,1);
          P_f_indirect=zeros(kz+2,1);
          
          P_f_direct(tl+stl+tempk,1)=1/n*trace(Sti);   P_f_direct(tl+stl+tempk+1,1)=1/n*trace(Sti*Wt*Sti*bhat(tempk));
          P_f_total(tl+stl+tempk,1)=1/n*sum(sum(Sti)); P_f_total(tl+stl+tempk+1,1)=1/n*sum(sum(Sti*Wt*Sti*bhat(tempk)));
          P_f_indirect(tl+stl+tempk,1)=P_f_total(tl+stl+tempk,1)-P_f_direct(tl+stl+tempk,1);
          
          sd_f_direct(tempk,tempt)=1/nt*P_f_direct'*SIGi1*P_f_direct;sd_f_direct(tempk,tempt)=sqrt(sd_f_direct(tempk,tempt));
          sd_f_total(tempk,tempt)=1/nt*P_f_total'*SIGi1*P_f_total;sd_f_total(tempk,tempt)=sqrt(sd_f_total(tempk,tempt));
          sd_f_indirect(tempk,tempt)=1/nt*P_f_indirect'*SIGi1*P_f_indirect;sd_f_indirect(tempk,tempt)=sqrt(sd_f_indirect(tempk,tempt));
          
    end
end

results.f_direct=f_direct;
results.f_total=f_total;
results.f_indirect=f_indirect;

results.sd_f_direct=sd_f_direct;
results.sd_f_total=sd_f_total;
results.sd_f_indirect=sd_f_indirect;

