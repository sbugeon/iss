function [r_norm,x_ls,r] = get_spot_residual(o,A,b,code_indices)

%% Orthogonal Matching Pursuit (OMP)
% https://github.com/seunghwanyoo/omp
% Input for data of dimension n:
%   A: dictionary (matrix) [n x nAtoms]
%   b: signals [n x nData]
%   SpotWeightFactor(o.nRounds x nData): The SpotWeighting used to normalise influence of
%   rounds. Should be the value after background removal for all iterations. 
% Output:
%   r_norm: norm of residual after removing codes specified by
%   code_indices.
%   x_ls(:,i): coef vector for code A(:,code_indices(i)) [nData x nAtoms]

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

