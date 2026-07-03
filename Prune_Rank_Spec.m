function [U_r, SV_r] = Prune_Rank_Spec(B, tau)
% Prune_Rank_Spec  Spectral-norm rank pruning of the final QB block (Algorithm 4.2).
%
%   [U_r, SV_r] = Prune_Rank_Spec(B, tau)
%
%   Given the final retained block B (b-by-n) and the target tolerance tau,
%   returns the truncated SVD factors that keep exactly the singular
%   directions of B whose singular value is at least tau:
%
%       ||B - U_r*SV_r||_2 = sigma_{r+1}(B) < tau.
%
%   The retained part is returned in factored form as U_r (left singular
%   vectors) and SV_r = S_r*V_r' (scaled right singular vectors), so that the
%   pruned block is simply B_r = U_r*SV_r.  Since S_r is diagonal, SV_r is
%   formed by scaling the rows of V_r' by the retained singular values sv(r)
%   -- no diagonal matrix is ever built.
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
%     SV_r  r-by-n scaled right factor S_r*V_r' (so B_r = U_r*SV_r)
%
%   The retained rank r satisfies sigma_r(B) >= tau > sigma_{r+1}(B).

    [U, S, V] = svd(B, 'econ');
    sv = diag(S);

    r = find(sv >= tau);        % singular values are sorted, so r is a leading block

    U_r  = U(:, r);
    SV_r = sv(r) .* V(:, r)';   % S_r*V_r' via row scaling (implicit expansion)
end
