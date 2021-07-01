function [time] = ndte_all(xx)
% function to compute NDTE for any of the eight conditions with surrogates for dbs80 
% call with number of rest or tasks to be computed
parcel='dbs80';
tasks={'REST1'; 'WM';'EMOTION'; 'SOCIAL'; 'MOTOR'; 'LANGUAGE'; 'GAMBLING'; 'RELATIONAL'};
direcs={'LR'};

% number of surrogate timeseries
ITER=100;

% open the local cluster profile
p = parcluster('local');

% open the parallel pool, recording the time it takes
time_pool = tic;
parpool(p,n);
Rtime_pool = toc(time_pool);
fprintf('Opening the parallel pool took %g seconds.\n', time_pool)

% output participant no
basein=['./ndte/' parcel '_data'];
baseout=[parcel '_processed'];
baseout='.';

% to make sure logdet is included in batch job
a=logdet(10);

% datafile
group=['hcp1003_' tasks{xx} '_LR_dbs80.mat'];
partic=[basein '/' group];

% load data for these participants
load(partic)
% number of brain regions in parcellation
N_areas=size(subject{1}.signal_filt,1);
Tmax=size(subject,2)

%% NDTE computation
for s=1:Tmax

  if (isnumeric(subject{s}))
    g{s}=nan;
  else
    XX=subject{s}.signal_filt;
    N=size(XX,1);

    try
        [GCval,GCr,Pval] = ndte_surrogates_fixlags_cs(XX,ITER);
    catch
        warning(['participant data does not exist: ' num2str(s)]);  
        [GCval,GCr,Pval] = nan;
    end

    g{s}.GCval=GCval;
    g{s}.GCr=GCr;
    g{s}.Pval=Pval;
    g{s}.NLags=10*ones(N,N);
  end
end;

% save result
save ([baseout '/GCall_' tasks{xx} '_ITER' num2str(ITER)],'g','-v7.3');

