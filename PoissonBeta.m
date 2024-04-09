function [samples_Pij, samples_Si, samples_Koni, samples_Koffi, Qiter] = ...
    PoissonBeta(count_data, size_factor, transcript_length, Tmax)
% PoissonBeta
% 
% Input: 
% 1. count_data: count data
% 2. size_factor
% 3. transcript_length
% 4. Tmax: a maximum number of Gibbs sampling interations
% 
% By Jong Kyoung Kim 
% jkkim@ebi.ac.uk
% Last Update: 10 July 2012

RandStream.setDefaultStream ...
     (RandStream('mt19937ar','seed',sum(100*clock)));

burnin = round(Tmax/2);

[num_gene, num_replicate] = size(count_data);
sum_count = sum(count_data, 2);

%% Initialization for hyperparameters
[C_param, A_param] = PoissonBetaInit(count_data, size_factor, transcript_length, num_gene, num_replicate);

%% Gibbs sampling
Qiter = zeros(Tmax,1);
samples_Pij = zeros(C_param.num_gene, C_param.num_replicate);
samples_Si = zeros(C_param.num_gene, Tmax - burnin);
samples_Koni = zeros(C_param.num_gene, Tmax - burnin);
samples_Koffi = zeros(C_param.num_gene, Tmax - burnin);
for t=1:Tmax
    fprintf('%s\n', 'Update Pij');
%     C_param = PoissonBetaUpdatePij(count_data, C_param, A_param);
    C_param.Pij = PoissonBetaUpdatePij_C(count_data, C_param, A_param);
    C_param.lnPij = log(C_param.Pij);
    C_param.lnoneminusPij = log(1-C_param.Pij);
    fprintf('%s\n', 'Update Si');
%     C_param = PoissonBetaUpdateSi(sum_count, A_param, C_param);
    C_param.Si = PoissonBetaUpdateSi_C(sum_count, A_param, C_param);
    C_param.lnSi = log(C_param.Si);
    fprintf('%s\n', 'Update Koni');
%     C_param = PoissonBetaUpdateKoni(A_param, C_param);
    C_param.Koni = PoissonBetaUpdateKoni_C(A_param, C_param);
    C_param.lnKoni = log(C_param.Koni);
    fprintf('%s\n', 'Update Koffi');
%     C_param = PoissonBetaUpdateKoffi(A_param, C_param);   
    C_param.Koffi = PoissonBetaUpdateKoffi_C(A_param, C_param);  
    C_param.lnKoffi = log(C_param.Koffi);
    Qiter(t) = PoissonBetaLogPosterior(count_data, A_param, C_param);
    fprintf('%s\t%d\t%d\n', 'PoissonBeta Gibbs Sampling Iteration', t, Qiter(t));
    if t > burnin
        samples_Pij = samples_Pij + C_param.Pij;
        samples_Si(:, t-burnin) = C_param.Si;
        samples_Koni(:, t-burnin) = C_param.Koni;
        samples_Koffi(:, t-burnin) = C_param.Koffi;
    end;
end;
samples_Pij = samples_Pij/(Tmax-burnin);