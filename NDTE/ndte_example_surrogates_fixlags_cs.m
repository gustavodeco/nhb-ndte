function [GCval,GCr,Pvals] = ndte_example_surrogates_fixlags_cs(X,ITER)
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
        Iyxz1=Iy-(logdet(cov(z))-logdet(cov(z(:,2:end))));
        Iypast=Hy1+logdet(cov(z(:,2:end)))-logdet(cov(z));
        Iyxz=Iyxz1/Iypast;
        if Iyxz<=0
          Iyxz=1e-15;
        end
        
        for iter=1:ITER
            drawnbin = randi([ceil(Tmax*0.05) ceil(Tmax*0.95)]);
            auxbin = [drawnbin:Tmax 1:drawnbin-1];
            XSi = X(i,auxbin);
            drawnbin = randi([ceil(Tmax*0.05) ceil(Tmax*0.95)]);
            auxbin = [drawnbin:Tmax 1:drawnbin-1];
            XSk = X(k,auxbin);
            for j=MaxLag+1:Tmax
                y_presuri(j-MaxLag,1:NLags(i,i))=XSi(j:-1:j-NLags(i,i)+1);
                y_presur(j-MaxLag,1:NLags(i,k))=XSk(j:-1:j-NLags(i,k)+1);
            end
            ys=squeeze(y_presuri(:,1:NLags(i,i)));
            y1s=squeeze(y_presuri(:,1));
            Hy1s=log(var(y1s));
            Iys=(logdet(cov(ys))-logdet(cov(ys(:,2:end))));
            xs=y_presur(:,2:NLags(i,k));
            if size(xs,1)==1
                xs=xs';
            end
            zs=horzcat(ys,xs);
            Iyxz1s=Iys-(logdet(cov(zs))-logdet(cov(zs(:,2:end))));
            Iypasts=Hy1s+logdet(cov(zs(:,2:end)))-logdet(cov(zs));
            Iyxzs(iter)=Iyxz1s/Iypasts;
            if Iyxzs(iter)<=0
               Iyxzs(iter)=1e-15;
            end
        end
        GCr(i,k)=(Iyxz-mean(Iyxzs))/std(Iyxzs);
        maxrange=max(max(Iyxzs),Iyxz);
        range=0:maxrange/1000:maxrange;
        f=ksdensity(Iyxzs,range,'Support','positive','Function','cdf');
        [aux idx]=min(abs(range-Iyxz));
        Pvals(i,k)=1-f(idx);
        GCval(i,k)=Iyxz;
    end
end