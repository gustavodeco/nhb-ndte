function [GCsim] = pair_granger_norm(X,NLags)

N=size(X,1);
MaxLag=max(max(NLags));
Numdata=size(X,2)-MaxLag;
y_pre=zeros(N,Numdata,MaxLag+1);

for i=1:N
    for p=1:MaxLag
        for j=MaxLag+1:size(X,2)
            y_pre(i,j-MaxLag,1:p+1)=X(i,j:-1:j-p);
        end
    end
end

GCsim=zeros(N,N);

for i=1:N
    nodelist=1:N;
    nodelist(i)=[];
    y=squeeze(y_pre(i,:,1:NLags(i,i)));
    y1=squeeze(y_pre(i,:,1));
    Hy1=log(var(y1));
    Iy=(logdet(cov(y))-logdet(cov(y(:,2:end))));
    for k=nodelist
        x=squeeze(y_pre(k,:,2:NLags(i,k)));
        if size(x,1)==1
            x=x';
        end
        z=horzcat(y,x);
        Iyxz=Iy-(logdet(cov(z))-logdet(cov(z(:,2:end))));
        Iypast=Hy1+logdet(cov(z(:,2:end)))-logdet(cov(z));
        GCsim(i,k)=Iyxz/Iypast;
    end
end