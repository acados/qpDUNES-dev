
%% initialize

% instructions:
% - compile blasfeo with desired LA
% - set __USE_BLASFEO_FLAG__ in qpDUNES
% - make qpDUNES
% - run C code (results stored in txt files)
% - update LA field below {'BLASFEO_HP', 'BLASFEO_RF', 'BLASFEO_WR_NETLIB', 'BLASFEO_WR_OPENBLAS','QPDUNES'}
% - update NM (as in C code)
% - run script to convert timings to mat files

clear variables; close all; clc;

% choose maximum number of iterations (to allocate memory)
MAXIT = 20;

% choose experiments (number of masses)
NM = [3 4 5 6 7 8 9 10 11 12 13 14 15];

% choose linear algebra used in experiments:
LA = 'UNKNOWN';

%% read in timings

% initialize iteration timings 
tNwtnSetup  = NaN(MAXIT, numel(NM));
tNwtnFactor = NaN(MAXIT, numel(NM));
tNwtnSolve  = NaN(MAXIT, numel(NM));
tStageQps   = NaN(MAXIT, numel(NM));
tLineSearch = NaN(MAXIT, numel(NM));
tOverhead   = NaN(MAXIT, numel(NM));
tIter       = NaN(MAXIT, numel(NM));
nIter       = NaN(1, numel(NM));

% load from txt files
for ii = 1:length(NM)
    load(['tNwtnSetup_' num2str(NM(ii)) '.txt']);
    load(['tNwtnFactor_' num2str(NM(ii)) '.txt']);
    load(['tNwtnSolve_' num2str(NM(ii)) '.txt']);
    load(['tQP_' num2str(NM(ii)) '.txt']);
    load(['tLineSearch_' num2str(NM(ii)) '.txt']);
    load(['tExtra_' num2str(NM(ii)) '.txt']);
    load(['tIter_' num2str(NM(ii)) '.txt']);
    
    eval(['nIter(ii) = length(tNwtnSetup_' num2str(NM(ii)) ');' ]);    
    eval(['tNwtnSetup(1:nIter(ii), ii) = tNwtnSetup_'  num2str(NM(ii)) ';'])
    eval(['tNwtnFactor(1:nIter(ii), ii) = tNwtnFactor_'  num2str(NM(ii)) ';'])
    eval(['tNwtnSolve(1:nIter(ii), ii) = tNwtnSolve_'  num2str(NM(ii)) ';'])
    eval(['tStageQps(1:nIter(ii), ii) = tQP_'  num2str(NM(ii)) ';'])
    eval(['tLineSearch(1:nIter(ii), ii) = tLineSearch_'  num2str(NM(ii)) ';'])
    eval(['tOverhead(1:nIter(ii), ii) = tExtra_'  num2str(NM(ii)) ';'])
    eval(['tIter(1:nIter(ii), ii) = tIter_'  num2str(NM(ii)) ';'])

end

%% sanity checks

percOverhead = 100*(tOverhead ./ tIter);

if any(percOverhead > 10)
    error('too much overhead detected in timings!')
end

if any(nIter > MAXIT)
    error('MAXIT is too small!');
end

plot(2*NM, 1000*(tNwtnFactor(1,:) + tNwtnSolve(1,:) ),'-o','linewidth',2)

%% save processed data

strngs = {'NwtnSetup', 'NwtnFactor', 'NwtnSolve', 'StageQps', 'LineSearch', 'Overhead', 'Iter' };

for ii = 1:length(strngs)

    eval(['t' strngs{ii} '_' LA ' = t' strngs{ii} ';']);
    eval(['save(''t' strngs{ii} '_' LA ''', ''t' strngs{ii} '_' LA ''');']);  
end
