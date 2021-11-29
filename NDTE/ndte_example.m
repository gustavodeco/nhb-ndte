% compute NDTE for any of the eight conditions with surrogates for dbs80 

% number of surrogate timeseries
ITER=100;

% to make sure logdet is included in batch job
a=logdet(10);

% datafile
% for testing purposes
group='hcpunrelated100_EMOTION.mat';
partic=[group];

% load data for these participants
load(partic)
% number of brain regions in parcellation
N_areas=size(subject{1}.dbs80ts,1);
N_subjects = size(subject,2);

%% NDTE computation
for s=1:N_subjects

    TR=0.72;  % Repetition Time (seconds)

    % Bandpass filter settings
    fnq=1/(2*TR);                 % Nyquist frequency
    flp = 0.008;                    % lowpass frequency of filter (Hz)
    fhi = 0.08;                    % highpass
    Wn=[flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency
    k=2;                          % 2nd order butterworth filter
    [bfilt,afilt]=butter(k,Wn);   % construct the filter

    BOLD = subject{s}.dbs80ts;
    
    for seed=1:80
        BOLD(seed,:)=zscore(detrend(BOLD(seed,:)));
        subject{s}.signal_filt(seed,:) = filtfilt(bfilt,afilt,BOLD(seed,:));
    end

    XX=subject{s}.signal_filt;

    N=size(XX,1);

    [GCval,GCr,Pval] = ndte_example_surrogates_fixlags_cs(XX,ITER);

    g{s}.GCval=GCval;
    g{s}.GCr=GCr;
    g{s}.Pval=Pval;
    g{s}.NLags=10*ones(N,N);
end


% save result
%save ([baseout '/GCall_' tasks{xx} '_ITER' num2str(ITER)],'g','-v7.3');

