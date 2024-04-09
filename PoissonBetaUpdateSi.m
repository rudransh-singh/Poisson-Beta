function C_param = PoissonBetaUpdateSi(sum_count, A_param, C_param)
% PoissonBetaUpdateSi
% 
% By Jong Kyoung Kim 
% jkkim@ebi.ac.uk
% Last Update: 02 July 2012

for ii=1:C_param.num_gene
    sum_Pij = sum(A_param.size_factor(ii,:).*(C_param.Pij(ii,:)-1));
    sum_x = sum_count(ii);
    f = @(x) -1*x./A_param.beta_si(ii) + (A_param.alpha_si(ii)-1).*log(x) + ...
        sum_x*log(x) + x*sum_Pij;
    C_param.Si(ii) = PoissonBetaSliceSampleGamma(f, C_param.Si(ii), C_param.Si(ii)/2);
end;
C_param(1).lnSi = log(C_param(1).Si);