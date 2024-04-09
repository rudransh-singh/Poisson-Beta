function C_param = PoissonBetaUpdateKoni(A_param, C_param)
% PoissonBetaUpdateKoni
% 
% By Jong Kyoung Kim 
% jkkim@ebi.ac.uk
% Last Update: 02 July 2012

for ii=1:C_param.num_gene
    f = @(x) -1*x./A_param.beta_koni(ii) + (A_param.alpha_koni(ii)-1).*log(x) + ...
        sum(gammaln(x + C_param.Koffi(ii)) - gammaln(x) + (x-1)*C_param.lnoneminusPij(ii,:));
    C_param.Koni(ii) = PoissonBetaSliceSampleGamma(f, C_param.Koni(ii), C_param.Koni(ii)/2);
end;
C_param(1).lnKoni = log(C_param(1).Koni);