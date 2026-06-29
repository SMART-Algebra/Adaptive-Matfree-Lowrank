function s = Prune_Rank_Fro(B, alpha, tau)
% Prune_Rank_Fro  Frobenius-norm rank pruning of the final QB block (Algorithm 3.2).
%
%   s = Prune_Rank_Fro(B, alpha, tau)
%
%   Given the final block B (b-by-n), the current residual-norm estimate
%   alpha = ||E||_F, and the target tolerance tau, returns the index set s of
%   the rows/columns of the block to KEEP so that the pruned approximation
%   still satisfies the tolerance.
%
%   The selection uses the bottom-up accumulation
%       ||E||_F^2 + sum_{j in r} ||b_j||_2^2 < tau^2,
%   adding the smallest row contributions first.  This avoids the subtractive
%   cancellation of top-down (energy-subtraction) schemes.
%
%   Inputs
%     B      b-by-n final block (rows are b_j = e_j' * B)
%     alpha  current residual Frobenius-norm estimate (>= 0)
%     tau    target tolerance (absolute, i.e. tau = epsilon or epsilon*||A||_F)
%
%   Output
%     s      indices in {1,...,b} to retain (complement of the discarded set)

    b = size(B, 1);

    % Row norms of the block, sorted ascending so we discard the
    % lowest-energy rows first.
    w = vecnorm(B, 2, 2);
    [w_sorted, prm] = sort(w, 'ascend');

    r = [];                                 % indices to discard
    for j = 1:b
        alpha_new = sqrt(alpha^2 + w_sorted(j)^2);
        if alpha_new >= tau
            break                           % discarding this row would break tau
        end
        alpha = alpha_new;
        r = [r, prm(j)];                    %#ok<AGROW>
    end

    s = setdiff(1:b, r, 'stable');          % keep everything not discarded
end
