function C_param = PoissonBetaUpdatePij(count_data, C_param, A_param)
% PoissonBetaUpdatePij
% 
% By Jong Kyoung Kim 
% jkkim@ebi.ac.uk
% Last Update: 02 July 2012

for ii=1:C_param.num_gene
    for jj=1:C_param.num_replicate
        f = @(x) (C_param.Koffi(ii)-1).*log(x) + ...
            (C_param.Koni(ii)-1).*log(1-x) + ...
            log(1-x)*count_data(ii,jj) + ...
            A_param.size_factor(ii,jj)*(C_param.Si(ii)*(x-1));
        C_param.Pij(ii,jj) = PoissonBetaSliceSampleBeta(f, C_param.Pij(ii,jj), C_param.Pij(ii,jj)/2);
    end;
end;
C_param(1).lnPij = log(C_param(1).Pij);
C_param(1).lnoneminusPij = log(1-C_param(1).Pij);