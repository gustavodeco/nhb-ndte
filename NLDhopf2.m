function fval=NLDhopf2(x)
global N GCdata GCsim NLags we dt Tmax sig2 tau Iext bfilt2 afilt2 TR C Cnew MGC NSUB a omega lagh;
nn=0;
Cnew=zeros(N,N);
for i=1:N
    for j=1:N
        if (C(i,j)>0 || j==N-i+1) 
            nn=nn+1;
            Cnew(i,j)=x(nn);
        end
    end
end

dsig = sqrt(dt*sig2);
dsig=repmat(dsig,1,2);
wC = we*Cnew;
sumC = repmat(sum(wC,2),1,2); % for sum Cij*xj
for nsub=1:NSUB
    xs=zeros(Tmax,N);
    %number of iterations, 100 willk�hrlich, weil reicht in diesem Fall
    z = 0.1*ones(N,2); % --> x = z(:,1), y = z(:,2)
    nn=0;
    % discard first 3000 time steps
    zpast=repmat(z,1,lagh);
    for t=0:dt:1000
        suma = wC*zpast(:,1:2) - sumC.*zpast(:,1:2); % sum(Cij*xi) - sum(Cij)*xj
        zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
        
        z = z + dt*((a.*z + zz.*omega - z.*(z.*z+zz.*zz))/tau + suma) + dsig.*randn(N,2);
        zpast=horzcat(zpast(:,3:end),z);
    end
    % actual modeling (x=BOLD signal (Interpretation), y some other oscillation)
    for t=0:dt:((Tmax-1)*TR)
        suma = wC*zpast(:,1:2) - sumC.*zpast(:,1:2); % sum(Cij*xi) - sum(Cij)*xj
        zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
        z = z + dt*((a.*z + zz.*omega - z.*(z.*z+zz.*zz))/tau + suma) + dsig.*randn(N,2);
        zpast=horzcat(zpast(:,3:end),z);
        if abs(mod(t,TR))<0.01
            nn=nn+1;
            xs(nn,:)=z(:,1)';
        end
    end
    %%%%
    BOLD=xs';
    signal_filt=zeros(N,nn);

    for seed=1:N
        BOLD(seed,:)=demean(detrend(BOLD(seed,:)));
        signal_filt(seed,:) =filtfilt(bfilt2,afilt2,BOLD(seed,:));
    end
    [GCsim] = pair_granger_norm(signal_filt,NLags);
    GCsim2(nsub,:,:)=GCsim;
end
GCsim=squeeze(mean(GCsim2));
xco=1-corrcoef(GCdata(:),GCsim(:));
fval=xco(2);
