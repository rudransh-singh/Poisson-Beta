function C_param = PoissonBetaUpdateKoffi(A_param, C_param)
% PoissonBetaUpdateKonc
% 
% By Jong Kyoung Kim 
% jkkim@ebi.ac.uk
% Last Update: 02 July 2012

for ii=1:C_param.num_gene
    f = @(x) -1*x./A_param.beta_koffi(ii) + (A_param.alpha_koffi(ii)-1).*log(x) + ...
        sum(gammaln(x + C_param.Koni(ii)) - gammaln(x) + ...
        (x-1)*C_param.lnPij(ii,:));
    C_param.Koffi(ii) = PoissonBetaSliceSampleGamma(f, C_param.Koffi(ii), C_param.Koffi(ii)/2);
end;
C_param(1).lnKoffi = log(C_param(1).Koffi);