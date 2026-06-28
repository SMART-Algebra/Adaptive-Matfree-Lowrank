function [Q, B, k, Err_ind, Exact_err, S] = randQB_MF_Fro(A, b, epsilon, compute_exact)
% randQB_MF_Fro  Matrix-free, adaptive low-rank QB approximation in the
%                Frobenius norm (Algorithm 3.1).
%
%   [Q, B, k, Err_ind, Exact_err, S] = randQB_MF_Fro(A, b, epsilon)
%   [...] = randQB_MF_Fro(A, b, epsilon, compute_exact)
%
%   Builds Q (orthonormal columns) and B = Q'*A block by block until the
%   randomized indicator certifies  ||A - Q*B||_F <= epsilon*||A||_F, then
%   prunes the final block to a near-minimal rank.  A is accessed only
%   through matrix-vector products (the stopping rule uses the relative
%   tolerance, i.e. tau = epsilon*||A||_F, estimated on the fly via S).
%
%   Inputs
%     A             m-by-n matrix
%     b             block size
%     epsilon       relative tolerance (0 < epsilon < 1)
%     compute_exact (optional, default false) if true, also return the exact
%                   residual norm at each iteration.  This forms the dense
%                   residual A - Q*B and is NOT matrix-free; use only for
%                   diagnostics/experiments.
%
%   Outputs
%     Q, B          QB factors with ||A - Q*B||_F <~ epsilon*||A||_F
%     k             computed rank, size(Q,2)
%     Err_ind       randomized error indicator per iteration
%     Exact_err     exact residual norm per iteration ([] unless compute_exact)
%     S             running estimate of ||A||_F per iteration

    if nargin < 4 || isempty(compute_exact)
        compute_exact = false;
    end

    [m, n]  = size(A);
    min_mn  = min(m, n);
    maxit   = ceil(min_mn / b);

    % Preallocate factors and per-iteration diagnostics.
    Q = zeros(m, maxit * b);
    B = zeros(maxit * b, n);
    ncol = 0;                       % number of columns currently stored in Q

    Err_ind = zeros(maxit, 1);
    S       = zeros(maxit, 1);
    if compute_exact
        Exact_err = zeros(maxit, 1);
    else
        Exact_err = [];
    end

    sum_Z = 0;                      % running sqrt(sum_k ||Z_k||_F^2)
    k = min_mn;                     % default rank if the tolerance is never met

    for i = 1:maxit

        % Random sketch of A.  Folding 1/sqrt(b) into Omega makes ||Y||_F a
        % direct (unbiased) estimate of the residual Frobenius norm.
        Omega = (1/sqrt(b)) * randn(n, b);
        Z = A * Omega;

        % Running estimate S(i) ~ ||A||_F for the relative stopping rule.
        sum_Z = hypot(sum_Z, norm(Z, 'fro'));
        S(i)  = sum_Z / sqrt(i);

        % Implicitly sketch the residual:  Y = (A - Q*B)*Omega.
        if ncol == 0
            Y = Z;
        else
            Qc = Q(:, 1:ncol);
            Bc = B(1:ncol, :);
            Y  = Z - Qc * (Bc * Omega);
        end

        if compute_exact
            if ncol == 0
                Exact_err(i) = norm(A, 'fro');
            else
                Exact_err(i) = norm(A - Qc * Bc, 'fro');
            end
        end

        % Randomized stopping criterion (relative tolerance).
        Err_ind(i) = norm(Y, 'fro');
        if Err_ind(i) < epsilon * S(i)
            % ---- Rank pruning of the final block (Prune_Rank_Fro) ----
            alpha = Err_ind(i);
            tau   = epsilon * S(i);

            B_end = B(ncol-b+1:ncol, :);        % rows of the last block B_{i-1}
            w = vecnorm(B_end, 2, 2);
            [w_sorted, prm] = sort(w, 'ascend');

            r = [];
            for j = 1:b
                alpha_new = sqrt(alpha^2 + w_sorted(j)^2);
                if alpha_new >= tau
                    break
                end
                alpha = alpha_new;
                r = [r, prm(j)];                %#ok<AGROW>
            end
            s = setdiff(1:b, r, 'stable');      % indices of the last block to keep

            keep = ncol - b;                    % columns kept from earlier blocks
            Q = [Q(:, 1:keep), Q(:, keep+s)];
            B = [B(1:keep, :);  B(keep+s, :)];

            Err_ind = Err_ind(1:i);
            S       = S(1:i);
            if compute_exact
                Exact_err = Exact_err(1:i);
            end
            k = size(Q, 2);
            return
        end

        % Orthonormalize the new block, then re-orthogonalize against Q.
        [Qi, ~] = qr(Y, 'econ');
        if ncol > 0
            Qc = Q(:, 1:ncol);
            [Qi, ~] = qr(Qi - Qc * (Qc' * Qi), 'econ');
        end
        Bi = Qi' * A;

        Q(:, ncol+1:ncol+b) = Qi;
        B(ncol+1:ncol+b, :) = Bi;
        ncol = ncol + b;
    end

    % Tolerance never met: return the full factorization.
    Q = Q(:, 1:ncol);
    B = B(1:ncol, :);
    k = min(ncol, min_mn);
end
