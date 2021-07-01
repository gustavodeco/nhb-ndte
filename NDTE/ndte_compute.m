function [Pvals GCval]=ndte_compute(X,y_pre,i,ITER,MaxLag)
    
    N=size(X,1);
    Tmax=size(X,2);
    NLags=10*ones(N,N);

	Pvals=zeros(1,N);
    GCval=zeros(1,N);

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
        
        % create surrogate timeseries
        for iter=1:ITER
            Iyxzs(iter)=ndte_surrogates_fixlags_cs(X,i,k,MaxLag,iter);
        end
        
        GCr(i,k)=(Iyxz-mean(Iyxzs))/std(Iyxzs);
        maxrange=max(max(Iyxzs),Iyxz);
        range=0:maxrange/1000:maxrange;
        f=ksdensity(Iyxzs,range,'Support','positive','Function','cdf');
        [aux idx]=min(abs(range-Iyxz));
        Pvals(k)=1-f(idx);
        GCval(k)=Iyxz;
    end
