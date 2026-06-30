function [U_r, S_r, V_r] = Prune_Rank_Spec(B, tau)
% Prune_Rank_Spec  Spectral-norm rank pruning of the final QB block (Algorithm 4.2).
%
%   [U_r, S_r, V_r] = Prune_Rank_Spec(B, tau)
%
%   Given the final retained block B (b-by-n) and the target tolerance tau,
%   returns the truncated SVD factors that keep exactly the singular
%   directions of B whose singular value is at least tau:
%
%       ||B - U_r*S_r*V_r'||_2 = sigma_{r+1}(B) < tau.
%
%   Unlike the Frobenius-norm pruning (Prune_Rank_Fro), which discards whole
%   rows by their norm, the spectral-norm truncation error is exactly the
%   first discarded singular value (Eckart-Young), so we truncate B by its
%   singular value distribution.  Since B is short-and-wide (b small), its
%   SVD costs only O(n b^2).
%
%   Inputs
%     B     b-by-n final block
%     tau   target tolerance (absolute, i.e. tau = epsilon or epsilon*||A||_2)
%
%   Outputs
%     U_r   b-by-r retained left singular vectors
%     S_r   r-by-r diagonal matrix of retained singular values
%     V_r   n-by-r retained right singular vectors
%
%   The retained rank r satisfies sigma_r(B) >= tau > sigma_{r+1}(B).

    [U, S, V] = svd(B, 'econ');
    sv = diag(S);

    r = find(sv >= tau);        % singular values are sorted, so r is a leading block

    U_r = U(:, r);
    S_r = S(r, r);
    V_r = V(:, r);
end
