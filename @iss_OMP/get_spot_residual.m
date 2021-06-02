function [r_norm,x_ls,r] = get_spot_residual(o,A,b,code_indices)

% Orthogonal Matching Pursuit (OMP)
% https://github.com/seunghwanyoo/omp
% Input for data of dimension n:
%   A: dictionary (matrix) [n x nAtoms]
%   b: signals [n x nData]
%   code_indices: atoms at these indices will be found.
% Output:
%   r_norm: norm of residual after removing codes specified by
%   code_indices.
%   x_ls(:,i): coef vector for code A(:,code_indices(i)) [nData x nAtoms]
A_omega = A(:,code_indices);
x_ls = A_omega \ b;  % Aomega * x_ls = b
r = b - A_omega * x_ls; % get residual
r_norm = vecnorm(r,2,1)';
x_ls = x_ls';

end

