function [U, S, V, sv_drop] = QB_to_SVD(Q, B, trunc_type, trunc_val)
% QB_to_SVD  Convert a QB approximation A ~ Q*B into a (truncated) SVD.
%
%   [U, S, V]          = QB_to_SVD(Q, B)
%   [U, S, V, sv_drop] = QB_to_SVD(Q, B, trunc_type, trunc_val)
%
%   Given Q (m-by-k, orthonormal columns) and B (k-by-n) from randQB_MF_Fro
%   or randQB_MF_Spec, returns U (m-by-r), S (r-by-r diagonal) and V (n-by-r)
%   with
%
%       Q*B  ~  U*S*V',      U'*U = I,   V'*V = I,   r <= k.
%
%   The factors come from the SVD of the small matrix B = U_B*S*V' and
%   U = Q*U_B.  The cost is O((m+n) k^2) flops and NO products with A: this
%   is a post-processing step.  Without truncation (default) r = k and
%   U*S*V' = Q*B up to rounding.
%
%   Optional truncation of the trailing singular values s_j of B:
%     trunc_type   'none'  (default) keep all k triplets
%                  'rank'  keep the leading r = trunc_val triplets
%                  'fro'   drop trailing s_j while  sqrt(sum of dropped s_j^2) <= trunc_val
%                  'spec'  drop the trailing s_j with  s_j <= trunc_val
%     trunc_val    the rank (for 'rank') or an ABSOLUTE threshold (for 'fro'
%                  and 'spec').  Since ||B||_F = ||Q*B||_F and ||B||_2 = ||Q*B||_2,
%                  a relative threshold is trunc_val = rel*norm(B,'fro') or
%                  rel*norm(B).
%
%   Outputs
%     U, S, V      factors of the (truncated) approximation; s = diag(S)
%     sv_drop      the dropped singular values of B (empty if none), so that
%                      ||Q*B - U*S*V'||_F = norm(sv_drop),
%                      ||Q*B - U*S*V'||_2 = max(sv_drop).
%
%   Error of the truncated result.  The randQB_MF_* routines return B = Q'*A,
%   so A - Q*B is orthogonal to range(Q) and the two error terms add in
%   quadrature:
%       ||A - U*S*V'||_F^2  =  ||A - Q*B||_F^2 + norm(sv_drop)^2,
%       ||A - U*S*V'||_2   <=  sqrt( ||A - Q*B||_2^2 + max(sv_drop)^2 ).
%   The Frobenius-norm routines prune the final block row by row, so a 'fro'
%   truncation at the leftover budget can lower the rank further at no cost.
%
%   Pseudoinverse.  U and V have orthonormal columns and S is diagonal, so
%       pinv(U*S*V') = V * (S \ U')                       % r-by-r diagonal solve
%   and the Tikhonov-regularized version with parameter alpha is
%       V * ((S.^2 + alpha*eye(r)) \ (S*U')),
%   neither of which forms an m-by-n matrix.
%
%   Complex Q and B are handled (' is the conjugate transpose throughout).

    if nargin < 3 || isempty(trunc_type)
        trunc_type = 'none';
    end

    k = size(Q, 2);
    if size(B, 1) ~= k
        error('QB_to_SVD:dims', ...
              'size(B,1) = %d must equal size(Q,2) = %d.', size(B, 1), k);
    end

    % SVD of the small factor: B = U_B * S * V'.
    [U_B, S, V] = svd(B, 'econ');
    sv = diag(S);                       % sorted in descending order
    kmin = numel(sv);

    switch lower(trunc_type)
        case 'none'
            r = kmin;
        case 'rank'
            r = min(max(round(trunc_val), 0), kmin);
        case 'fro'
            % Drop the smallest singular values first, accumulating bottom-up
            % (the same accumulation as Prune_Rank_Fro: no subtractive
            % cancellation).  tail(j) = norm of the j smallest values.
            tail  = sqrt(cumsum(flipud(sv).^2));
            ndrop = find(tail <= trunc_val, 1, 'last');
            if isempty(ndrop)
                ndrop = 0;
            end
            r = kmin - ndrop;
        case 'spec'
            r = nnz(sv > trunc_val);    % sorted, so the kept set is a leading block
        otherwise
            error('QB_to_SVD:trunc_type', 'Unknown trunc_type ''%s''.', trunc_type);
    end

    sv_drop = sv(r+1:end);
    U = Q * U_B(:, 1:r);
    S = diag(sv(1:r));
    V = V(:, 1:r);
end
