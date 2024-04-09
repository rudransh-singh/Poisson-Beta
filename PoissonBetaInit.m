function [C_param, A_param] = PoissonBetaInit(count_data, size_factor, transcript_length, num_gene, num_replicate)
% PoissonBetaInit
% 
% Input:
% 1. count_data
% 2. size_factor
% 3. transcript_length
% 4. num_gene: a number of genes
% 5. num_replicate: a number of replicates
% 
% By Jong Kyoung Kim 
% jkkim@ebi.ac.uk
% Last Update: 25 September 2012

A_param = struct('alpha_si', {}, 'beta_si', {}, 'alpha_koni', {}, ...
    'beta_koni', {}, 'alpha_koffi', {}, 'beta_koffi', {}, ...
    'count_data_gammaln', {}, 'count_data_row_gammaln', {}, ...
    'size_factor', {}, 'ln_size_factor', {}, 'ln_size_factor_sum', {}, 'count_ln_size_factor', {});
A_param(1).size_factor = repmat(size_factor, num_gene, 1).*repmat(transcript_length, 1, num_replicate);
A_param(1).alpha_si = ones(num_gene,1);
A_param(1).beta_si = max(count_data./A_param.size_factor,[],2);       
A_param(1).alpha_koni = ones(num_gene,1); 
A_param(1).beta_koni = 100*ones(num_gene,1); 
A_param(1).alpha_koffi = ones(num_gene,1); 
A_param(1).beta_koffi = 100*ones(num_gene,1); 
A_param(1).count_data_gammaln = sum(sum(gammaln(count_data+1)));
A_param(1).count_data_row_gammaln = zeros(num_gene,1);
A_param(1).ln_size_factor = log(A_param(1).size_factor);
A_param(1).ln_size_factor_sum = sum(sum(count_data.*A_param.ln_size_factor));
for ii=1:num_gene
    A_param(1).count_data_row_gammaln(ii) = sum(gammaln(count_data(ii,:)+1));
    A_param(1).count_ln_size_factor(ii) = sum(count_data(ii,:).*A_param(1).ln_size_factor(ii,:));
end;
%% Initialization for variational distributions

C_param = struct('num_gene', {}, 'num_replicate', {}, 'Pij', {}, 'Si', {}, ...
    'Koni', {}, 'Koffi', {}, 'lnSi', {}, 'lnoneminusPij', {}, 'lnPij', {}, 'lnKoni', {}, ...
    'lnKoffi', {});
C_param(1).num_gene = num_gene;
C_param(1).num_replicate = num_replicate;
C_param(1).Si = mean(mean(count_data./A_param.size_factor))*ones(num_gene,1);
C_param(1).Koni = rand(num_gene, 1);
C_param(1).Koffi = rand(num_gene, 1);
mean_size_factor = mean(size_factor);
for ii=1:num_gene
    x = count_data(ii, :);
    e1 = mean(x);
    e2 = sum(x.*(x-1))/length(x);
    e3 = sum(x.*(x-1).*(x-2))/length(x);    
    r1 = e1; r2 = e2/e1; r3 = e3/e2;
    kon_hat = 2*r1*(r3-r2)/(r1*r2-2*r1*r3+r2*r3);
    koff_hat = 2*(r2-r1)*(r1-r3)*(r3-r2)/((r1*r2-2*r1*r3+r2*r3)*(r1-2*r2+r3));
    s_hat = (-r1*r2+2*r1*r3-r2*r3)/(r1-2*r2+r3);
    s_hat = s_hat/(transcript_length(1)*mean_size_factor);
    if s_hat > 0 && kon_hat > 0 && koff_hat > 0
        C_param(1).Si(ii) = s_hat;
        C_param(1).Koni(ii) = kon_hat;
        C_param(1).Koffi(ii) = koff_hat;
    end;
end;
C_param(1).lnSi = log(C_param(1).Si);
C_param(1).lnKoni = log(C_param(1).Koni);
C_param(1).lnKoffi = log(C_param(1).Koffi);

C_param(1).Pij = 0.5*ones(num_gene, num_replicate);
C_param(1).lnPij = log(C_param(1).Pij);
C_param(1).lnoneminusPij = log(1-C_param(1).Pij);



