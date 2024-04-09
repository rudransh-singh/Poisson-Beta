function [pred_dist, xij] = PoissonBetaPostPredictive(samples_Pcj, samples_Sc, size_factor, cluster_index, x_min, x_max)
% PoissonBetaPostPredictive
% 
% By Jong Kyoung Kim 
% jkkim@ebi.ac.uk
% Last Update: September 07 2011

[num_gene, num_replicate] = size(size_factor);
[~, ~, T] = size(samples_Pcj);
xij = linspace(x_min, x_max, 10000);
pred_dist = zeros(length(xij),1);
for ii=1:num_gene
    for tt=1:T
        for jj=1:num_replicate
            f = @(x) exp(x*(log(size_factor(ii,jj)) + log(samples_Sc(cluster_index,tt)) + ...
                log(1-samples_Pcj(cluster_index, jj, tt))) + ...
                size_factor(ii,jj)*samples_Sc(cluster_index,tt)*(samples_Pcj(cluster_index, jj, tt)-1)-gammaln(x+1));
            pred_dist = pred_dist + f(xij');
        end;
    end;
end;




