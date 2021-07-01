function [GCval,GCr,Pvals] = ndte_surrogates_fixlags_cs(X, ITER)
MaxLag=10;
N=size(X,1);
Tmax=size(X,2);
NLags=10*ones(N,N);

Numdata=Tmax-MaxLag;
y_pre=zeros(N,Numdata,MaxLag+1);

for i=1:N
    for p=1:MaxLag
        for j=MaxLag+1:Tmax
            y_pre(i,j-MaxLag,1:p+1)=X(i,j:-1:j-p);
        end
    end
end

GCr=zeros(N,N);

parfor i=1:N
    [Pvals(i,:) GCval(i,:)]=ndte_compute(X,y_pre,i,ITER,MaxLag);    
end
