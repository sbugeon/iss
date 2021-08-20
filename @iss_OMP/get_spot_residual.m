function [r_norm,x_ls,r] = get_spot_residual(o,A,b,code_indices)
%% [r_norm,x_ls,r] = get_spot_residual(o,A,b,code_indices)
% Extended from: https://github.com/seunghwanyoo/omp
% This is run from call_spots_omp_initial when know which genes to remove
% and want to remove the same genes from every SpotColor.
%
% Input for data of dimension n = o.nBP*o.nRounds
%   o: iss object
%   A: dictionary of gene codes [n x nGenes]
%   b: signals i.e. NormSpotColors [n x nData]
%   code_indices: genes to be removed from signals (same for all signals).
%       [1 x nGenesToRemove]
% Output
%   r_norm: norm of residual after removing genes specified by
%       code_indices. [n x nData]
%   x_ls(:,i): coef vector for code A(:,code_indices(i)) [nData x nGenes]
%   r: residual after removing genes specified by code_indices [n x 1]
%
%% 
A_omega = A(:,code_indices);
%x_ls = sum(A_weight.*bWeight)./sum(WeightFactor.^2);
if length(code_indices)==1 && size(A_omega,2)==1
    %Do all pixels at same time if all same gene
    %LS solution: x_ls = (A_weight'*A_weight)^-1*A_weight'*bWeight
    %x_ls = diag(A_weight'*b_weight./sum(A_weight.^2,1))'; %TOO MUCH MEMORY
    x_ls = sum((A_omega.*b))./sum(A_omega.^2,1);
else
    x_ls = A_omega \ b;  % Aomega * x_ls = b
end
r = b - A_omega * x_ls; % get residual
r_norm = vecnorm(r,2,1)';
x_ls = x_ls';


end

