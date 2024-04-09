function Q = PoissonBetaLogPosterior(count_data, A_param, C_param)
% PoissonBetaLogPosterior
% 
% Input: 
% 1. sum_count: sum(count_data,2)
% 2. A_param: hyperparameters
% 3. C_param: variational parameters
% 
% By Jong Kyoung Kim 
% jkkim@ebi.ac.uk
% Last Update: 02 July 2012

% log p(X|Z,p, s)
Q = sum(sum(count_data.*(repmat(C_param.lnSi, 1, C_param.num_replicate) + ...
    C_param.lnoneminusPij))) + A_param(1).ln_size_factor_sum + ...
    sum(sum(A_param.size_factor.*(repmat(C_param.Si, 1, C_param.num_replicate).*(C_param.Pij-1)))) - ...
    A_param.count_data_gammaln;

% log p(p|k_off, k_on)
Q = Q + sum(sum(gammaln(repmat(C_param.Koni + C_param.Koffi, 1, C_param.num_replicate)) - ...
    gammaln(repmat(C_param.Koni, 1, C_param.num_replicate)) - gammaln(repmat(C_param.Koffi, 1, C_param.num_replicate)) + ...
    repmat(C_param.Koffi-1, 1, C_param.num_replicate).*C_param.lnPij + ...
    repmat(C_param.Koni-1, 1, C_param.num_replicate).*C_param.lnoneminusPij));

% log p(s)
Q = Q + sum(-1*C_param.Si./A_param.beta_si + (A_param.alpha_si-1).*C_param.lnSi - ...
    A_param.alpha_si.*log(A_param.beta_si) - gammaln(A_param.alpha_si));

% log p(kon)
Q = Q + sum(-1*C_param.Koni./A_param.beta_koni + (A_param.alpha_koni-1).*C_param.lnKoni - ...
    A_param.alpha_koni.*log(A_param.beta_koni) - gammaln(A_param.alpha_koni));

% log p(koff)
Q = Q + sum(-1*C_param.Koffi./A_param.beta_koffi + (A_param.alpha_koffi-1).*C_param.lnKoffi - ...
    A_param.alpha_koffi.*log(A_param.beta_koffi) - gammaln(A_param.alpha_koffi));