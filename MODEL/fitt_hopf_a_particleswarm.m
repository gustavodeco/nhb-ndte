clear all;
global N GCdata FCdata GCsim NLags we dt Tmax sig2 tau Iext bfilt2 afilt2 TR C Cnew MGC NSUB a omega lagh Isubdiag;
load  empiricalHCPrest.mat;

GCdata=GC;
NLags=Lags;

load SC_dbs80_dilatecortex.mat;
C=sc_dilated;
C=C/max(max(C))*0.2;

load hcpunrelated100_REST_dbs80.mat;
Tmax=1200;
N_areas=80;
N=80;
NSUB=100;

Isubdiag = find(tril(ones(N_areas),-1));

% Parameters of the data
TR=0.72;  % Repetition Time (seconds)

% Bandpass filter settings
fnq=1/(2*TR);                 % Nyquist frequency
flp = 0.01;                    % lowpass frequency of filter (Hz)
fhi = 0.2;                    % highpass
Wn=[flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency
k=2;                          % 2nd order butterworth filter
[bfilt2,afilt2]=butter(k,Wn);   % construct the filter
clear fnq flp fhi Wn k 

%%%%%%%%%%%%%%%%%%

for nsub=1:NSUB
    clear signal_filt_data;
    signaldata = subject{nsub}.dbs80ts;
    for seed=1:N
        signaldata(seed,:)=demean(detrend(signaldata(seed,:)));
        signal_filt_data(seed,:) =filtfilt(bfilt2,afilt2,signaldata(seed,:));
    end
    
    for i=1:N
        for j=1:N
            FCdata2(nsub,i,j)=corr2(signal_filt_data(i,:)',signal_filt_data(j,:)');
        end
    end
end
FCdata=squeeze(mean(FCdata2));

%%%%%%%%%%%%%%

insub=1;

for nsub=1:NSUB
    clear PowSpect PowSpect2;
    [N, Tmax]=size(subject{nsub}.dbs80ts);
    Isubdiag = find(tril(ones(N),-1));
    TT=Tmax;
    Ts = TT*TR;
    freq = (0:TT/2-1)/Ts;
    signaldata = subject{nsub}.dbs80ts;
    signaldata=signaldata(:,1:Tmax);

    nfreqs=length(freq);
    
    for seed=1:N
        x=detrend(demean(signaldata(seed,:)));
        ts =zscore(filtfilt(bfilt2,afilt2,x));
        pw = abs(fft(ts));
        PowSpect(:,seed,insub) = pw(1:floor(TT/2)).^2/(TT/TR);
    end
    insub=insub+1;
end

Power_Areas=mean(PowSpect,3);
for seed=1:N
    Power_Areas(:,seed)=gaussfilt(freq,Power_Areas(:,seed)',0.01);
end

[maxpowdata,index]=max(Power_Areas);
f_diff = freq(index);

clear PowSpect Power_Areas ;


%%%%%%%%%%%%%%%%%%%%%%%%

dt=0.1*TR/2;
a=-0.02*ones(N,2);
sig2=0.02*ones(N,2);
tau=1;
we=0.2;
omega = repmat(2*pi*f_diff',1,2); omega(:,1) = -omega(:,1);

Cnew=C;
lagh=1;


display('Optimized Particle Swarm we');
nvars1=1;
xinit1=we;
lb1=0.01;
ub1=4;
initpop1=0.02*randn(10,nvars1)+repmat(xinit1,10,1);
options = optimoptions('particleswarm', 'InitialSwarmMatrix',initpop1,'Display', 'iter','MaxIterations',10);
[x,fval] = particleswarm(@NLDhopf1,nvars1,lb1,ub1,options);
we=x;

for iter=1:3
    
    display('Optimized Particle Swarm a');
    
    nvars=N;
    xinit=a(:,1)';
    lb=-0.8*ones(1,nvars);
    ub=0.8*ones(1,nvars);
    initpop=0.02*randn(20,nvars)+repmat(xinit,20,1);
    options = optimoptions('particleswarm', 'InitialSwarmMatrix',initpop,'Display', 'iter','MaxIterations',10);
    [x,fval] = particleswarm(@NLDhopf2a,nvars,lb,ub,options);
    a=x'.*ones(N,2);
    
    display('Optimized Particle Swarm lag');
    
    for lag=1:20
        fitt0(lag)=NLDhopf0(lag);
    end
    [aux lagh]=min(fitt0)
    
    display('Optimized Particle Swarm sig2');
    nvars3=N;
    xinit3=sig2(:,1)';
    lb3=0.01*ones(1,N);
    ub3=0.2*ones(1,N);
    initpop3=0.02*randn(20,nvars3)+repmat(xinit3,20,1);
    options = optimoptions('particleswarm', 'InitialSwarmMatrix',initpop3,'Display', 'iter','MaxIterations',10);
    [x,fval] = particleswarm(@NLDhopf3,nvars3,lb3,ub3,options);
    sig2=x'
    
    display('Optimized Particle Swarm C');
    
    nvars=0;
    for i=1:N
        for j=1:N
            if (C(i,j)>0 || j==N-i+1)
                nvars=nvars+1;
                xinit(nvars)=Cnew(i,j);
            end
        end
    end
    lb=zeros(1,nvars);
    ub=0.2*ones(1,nvars);
    initpop=0.02*randn(20,nvars)+repmat(xinit,20,1);
    options = optimoptions('particleswarm', 'InitialSwarmMatrix',initpop,'Display', 'iter','MaxIterations',10);
    [x,fval] = particleswarm(@NLDhopf2,nvars,lb,ub,options);
    display('Optimized Particle Swarm');
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
    save Ceff_hopf.mat Cnew we a lagh sig2;
end
NLDhopf2(x)

