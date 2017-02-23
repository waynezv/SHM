function [D_out] = zscore_normalize2(D_in,mu_in,std_in)
[D, N] = size(D_in);
D_out = D_in-repmat(mu_in, [1 N]);
D_out = D_out./repmat(std_in, [1 N]);